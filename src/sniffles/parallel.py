#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 27.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

from dataclasses import dataclass
import gc
import math

import pysam

from sniffles import leadprov
from sniffles import cluster
from sniffles import sv
from sniffles import postprocessing
from sniffles import snf
from sniffles import util

@dataclass
class Task:
    """
    表示一个单独的task
    """
    id: int
    sv_id: int
    contig: str
    start: int
    end: int
    assigned_process_id: int=None
    lead_provider: object=None
    bam: object=None
    tandem_repeats: list=None
    genotype_svs: list=None

    def build_leadtab(self,config):
        """
        从bam中获取线索，构建leadtab。
        """
        assert(self.lead_provider==None)

        self.bam=pysam.AlignmentFile(config.input, config.input_mode, require_index=True) # 重新读取bam
        self.lead_provider=leadprov.LeadProvider(config,self.id*config.task_read_id_offset_mult) # 根据task id创建read id起始值
        externals=self.lead_provider.build_leadtab(self.contig,self.start,self.end,self.bam)
        return externals,self.lead_provider.read_count

    def call_candidates(self,keep_qc_fails,config):
        """
        迭代每一种sv，从leads计算cluster，再检测svcall。
        """
        candidates=[]
        # 迭代每种类型的sv
        for svtype in sv.TYPES:
            for svcluster in cluster.resolve(svtype,self.lead_provider,config,self.tandem_repeats):
                for svcall in sv.call_from(svcluster,config,keep_qc_fails,self):
                    if config.dev_trace_read!=False:
                        cluster_has_read=False
                        for ld in svcluster.leads:
                            if ld.read_qname==config.dev_trace_read:
                                cluster_has_read=True
                        if cluster_has_read:
                            import copy
                            svcall_copy=copy.deepcopy(svcall)
                            svcall_copy.postprocess=None
                            print(f"[DEV_TRACE_READ] [3/4] [Task.call_candidates] Read {config.dev_trace_read} -> Cluster {svcluster.id} -> preliminary SVCall {svcall_copy}")
                    candidates.append(svcall)

        self.coverage_average_fwd,self.coverage_average_rev=postprocessing.coverage(candidates,self.lead_provider,config) # 某个task区域的平均覆盖度
        self.coverage_average_total=self.coverage_average_fwd+self.coverage_average_rev
        return candidates

    def finalize_candidates(self,candidates,keep_qc_fails,config):
        """
        迭代每个svcall，执行一系列过滤，返回所有通过过滤的svcall。
        """
        passed=[]
        for svcall in candidates:
            svcall.qc=svcall.qc and postprocessing.qc_sv(svcall,config)
            if not keep_qc_fails and not svcall.qc:
                continue
            svcall.qc=svcall.qc and postprocessing.qc_sv_support(svcall,self.coverage_average_total,config)
            if not keep_qc_fails and not svcall.qc:
                continue

            postprocessing.annotate_sv(svcall,config)

            svcall.qc=svcall.qc and postprocessing.qc_sv_post_annotate(svcall,config)

            if config.dev_trace_read!=False:
                cluster_has_read=False
                for ld in svcall.postprocess.cluster.leads:
                    if ld.read_qname==config.dev_trace_read:
                        cluster_has_read=True
                if cluster_has_read:
                    import copy
                    svcall_copy=copy.deepcopy(svcall)
                    svcall_copy.postprocess=None
                    print(f"[DEV_TRACE_READ] [4/4] [Task.finalize_candidates] Read {config.dev_trace_read} -> Cluster {svcall.postprocess.cluster.id} -> finalized SVCall, QC={svcall_copy.qc}: {svcall_copy}")

            if not keep_qc_fails and not svcall.qc:
                continue

            # 清空svcall中保存的leads信息
            svcall.finalize() #Remove internal information (not written to output) before sending to mainthread for VCF writing
            passed.append(svcall)
        return passed

    def combine(self,config):
        """
        多样本snf文件联合检出
        """
        # 读取SNF文件头信息
        samples_headers_snf={}
        for snf_info in config.snf_input_info:
            snf_in=snf.SNFile(config,open(snf_info["filename"],"rb"),filename=snf_info["filename"])
            snf_in.read_header()
            samples_headers_snf[snf_info["internal_id"]]=(snf_info["filename"],snf_in.header,snf_in)

            if config.combine_close_handles:
                snf_in.close()

        svcalls=[] # 存储最后检出的svcall

        #block_groups_keep_threshold=5000
        #TODO: Parameterize
        bin_min_size=config.combine_min_size # 合并cand的窗口, 默认100bp.
        bin_max_candidates=max(25,int(len(config.snf_input_info)*0.5)) # SNF文件数目一半，最大25
        overlap_abs=config.combine_overlap_abs # 2500.

        # 读取样本内部ID
        sample_internal_ids=set()
        for sample_internal_id,(sample_filename,sample_header,sample_snf) in samples_headers_snf.items():
            sample_internal_ids.add(sample_internal_id)

        #
        # Load candidate SVs from all samples for each block separately and cluster them based on start position
        #
        # 迭代每个block, 获取所有样本的candidate svs，并按照位置聚类
        had=False
        candidates_processed=0
        svtypes_candidates_bins={svtype: {} for svtype in sv.TYPES}
        groups_keep={svtype:list() for svtype in sv.TYPES}

        # 迭代每个block
        for block_index in range(self.start,self.end+config.snf_block_size,config.snf_block_size): # 这里是因为start=0,即第一个block
            
            # 读取某个block下所有样本的数据
            # 读取所有样本在某个block的数据
            samples_blocks={}
            for sample_internal_id,(sample_filename,sample_header,sample_snf) in samples_headers_snf.items():
                blocks=sample_snf.read_blocks(self.contig,block_index) # 获取该block的数据
                samples_blocks[sample_internal_id]=blocks

            # 迭代每种sv
            for svtype in sv.TYPES:
                bins={} # 分bin存储所有样本某个sv类型的所有cand.
                #svcandidates=[]
                for sample_internal_id,(sample_filename,sample_header,sample_snf) in samples_headers_snf.items():
                    blocks=samples_blocks[sample_internal_id] # 某个block内某个样本的数据
                    if blocks==None:
                        continue
                    for block in blocks:
                        # block是字典对象, 包含了某个task在当前block内所有的cands.
                        for cand in block[svtype]:
                            #if config.combine_pass_only and (cand.qc==False or cand.filter!="PASS"):
                            #    continue

                            cand.sample_internal_id=sample_internal_id

                            # 按照新的bin大小存储cand.
                            bin=int(cand.pos/bin_min_size)*bin_min_size # 分100bp存储cand
                            if not bin in bins:
                                bins[bin]=[cand]
                            else:
                                bins[bin].append(cand)
                        candidates_processed+=len(block[svtype])
                
                # 某种类型的sv没有cand
                if len(bins)==0:
                    continue
                
                # 按条件给sv分组
                curr_bin=0
                size=0
                svcands=[]
                keep=groups_keep[svtype] # 初始为[].
                sorted_bins=sorted(bins)
                last_bin=sorted_bins[-1]
                for curr_bin in sorted_bins:
                    svcands.extend(bins[curr_bin])
                    size+=bin_min_size

                    # 满足条件时,触发sv分组
                    if (not config.combine_exhaustive and len(svcands) >= bin_max_candidates) or curr_bin == last_bin:
                        # 非详尽模式且svcands数目超过阈值, 或者到达最后bin
                        if len(svcands)==0:
                            # 到达最后bin但是没有svcands
                            size=0
                            continue
                        prevkept=len(keep)
                        svgroups=cluster.resolve_block_groups(svtype,svcands,keep,config) # 将svcands进行分组, 以sv.SVGroup表示
                        groups_call=[]
                        keep=[]
                        for group in svgroups:
                            # group是sv.SVGroup类型
                            # 记录不在group内的样本的覆盖度信息
                            coverage_bin=int(group.pos_mean/config.coverage_binsize_combine)*config.coverage_binsize_combine # 覆盖度bin
                            for non_included_sample in sample_internal_ids-group.included_samples:
                                if samples_blocks[non_included_sample]!=None and coverage_bin in samples_blocks[non_included_sample][0]["_COVERAGE"]:
                                    coverage=samples_blocks[non_included_sample][0]["_COVERAGE"][coverage_bin]
                                else:
                                    coverage=0
                                if non_included_sample in group.coverages_nonincluded:
                                    group.coverages_nonincluded[non_included_sample]=max(coverage,group.coverages_nonincluded[non_included_sample])
                                else:
                                    group.coverages_nonincluded[non_included_sample]=coverage

                            if abs(group.pos_mean - curr_bin) < max(size*0.5,overlap_abs):
                                # 当前group与当前bin距离很近
                                keep.append(group)
                            else:
                                groups_call.append(group)
                        svcalls.extend(sv.call_groups(groups_call,config,self))
                        size=0
                        svcands=[]

                groups_keep[svtype]=keep

        for svtype in groups_keep:
            svcalls.extend(sv.call_groups(groups_keep[svtype],config,self))

        return svcalls,candidates_processed

@dataclass
class Process:
    """
    表示每个进程。
    """
    id: int
    process: object=None
    pipe_main: object=None
    externals: list=None

def Main(proc_id,config,pipe):
    try:
        if config.dev_profile:
            import cProfile
            cProfile.runctx("Main_Internal(proc_id,config,pipe)",globals(),locals(),sort="cumulative")
        else:
            Main_Internal(proc_id,config,pipe)
    except Exception as e:
        pipe.send(["worker_exception",""])
        raise e

def Main_Internal(proc_id,config,pipe):
    """
    每个进程调用的目标函数。
    """
    tasks={}
    while True:
        command,arg=pipe.recv() # 从管道接收命令和task对象

        if command=="call_sample":
            # call_sample模式
            task=arg
            result={}

            if config.snf != None or config.no_qc:
                qc=False
            else:
                qc=True

            _,read_count=task.build_leadtab(config)
            svcandidates=task.call_candidates(qc,config)
            svcalls=task.finalize_candidates(svcandidates,not qc,config)

            # 是否输出所有的svcalls
            if config.no_qc:
                result["svcalls"]=svcalls
            else:
                result["svcalls"]=[s for s in svcalls if s.qc]

            # 输出snf文件
            if config.snf != None: # and len(svcandidates):
                # 每个task生成snf文件
                snf_filename=f"{config.snf}.tmp_{task.id}.snf"

                # 把候选的sv输出到snf文件中
                with open(snf_filename,"wb") as handle:
                    snf_out=snf.SNFile(config,handle)
                    for cand in svcandidates:
                        snf_out.store(cand)
                    snf_out.annotate_block_coverages(task.lead_provider) # 给block添加覆盖度信息
                    snf_out.write_and_index() # 把block数据序列化并输出
                    handle.close()
                
                result["snf_filename"]=snf_filename
                result["snf_index"]=snf_out.get_index()
                result["snf_total_length"]=snf_out.get_total_length()
                result["snf_candidate_count"]=len(svcandidates)
                result["has_snf"]=True
            else:
                # 不生成snf文件
                result["has_snf"]=False

            #if config.vcf != None:
            #    svcalls=task.finalize_candidates(svcandidates,config)
            #    result["svcalls"]=svcalls

            result["task_id"]=task.id
            result["processed_read_count"]=read_count
            result["coverage_average_total"]=task.coverage_average_total
            pipe.send(["return_call_sample",result]) # 通过管道把结果发送回主进程
            del task
            gc.collect() # 垃圾回收

        elif command=="genotype_vcf":
            # genotype_vcf模式
            task=arg
            result={}

            qc=False
            _,read_count=task.build_leadtab(config)
            svcandidates=task.call_candidates(qc,config=config)
            svcalls=task.finalize_candidates(svcandidates,not qc,config=config)

            # 按5kb窗口存储目标sv
            binsize=5000
            binedge=int(binsize/10)
            genotype_svs_svtypes_bins={svtype:{} for svtype in sv.TYPES} # 以svtype -> bin -> sv存储
            # 把目标svs分bin存储
            for genotype_sv in task.genotype_svs:
                # 迭代每个目标sv
                genotype_sv.genotype_match_sv=None
                genotype_sv.genotype_match_dist=math.inf

                if not genotype_sv.svtype in genotype_svs_svtypes_bins:
                    #TODO: Warn about unsupported SVTYPE
                    # 目标sv不在程序支持的sv类型中
                    continue
                
                # 考虑sv跨多个bin的可能性
                bins=[int(genotype_sv.pos/binsize)*binsize]
                if genotype_sv.pos%binsize < binedge:
                    bins.append((int(genotype_sv.pos/binsize)-1)*binsize)
                if genotype_sv.pos%binsize > binsize-binedge:
                    bins.append((int(genotype_sv.pos/binsize)+1)*binsize)

                for bin in bins:
                    if not bin in genotype_svs_svtypes_bins[genotype_sv.svtype]:
                         genotype_svs_svtypes_bins[genotype_sv.svtype][bin]=[]
                    genotype_svs_svtypes_bins[genotype_sv.svtype][bin].append(genotype_sv)

            # 迭代每个候选sv，判断是否与目标sv匹配
            for cand in svcandidates:
                bin=int(cand.pos/binsize)*binsize
                if not bin in genotype_svs_svtypes_bins[cand.svtype]:
                    continue
                if cand.svtype=="BND":
                    # BND类型
                    for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                        dist=abs(genotype_sv.pos - cand.pos)
                        #if minlen>0 and dist < genotype_sv.genotype_match_dist and dist <= config.cluster_merge_bnd * 2:
                        if dist < genotype_sv.genotype_match_dist and dist <= config.cluster_merge_bnd:
                            if cand.bnd_info.mate_contig==genotype_sv.bnd_info.mate_contig:
                                genotype_sv.genotype_match_sv=cand
                                genotype_sv.genotype_match_dist=dist
                else:
                    # 非BND类型
                    for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                        # 判断sv的位置和大小
                        dist=abs(genotype_sv.pos - cand.pos) + abs(abs(genotype_sv.svlen) - abs(cand.svlen))
                        minlen=float(min(abs(genotype_sv.svlen),abs(cand.svlen)))
                        if minlen>0 and dist < genotype_sv.genotype_match_dist and dist <= config.combine_match * math.sqrt(minlen) and dist <= config.combine_match_max:
                            genotype_sv.genotype_match_sv=cand
                            genotype_sv.genotype_match_dist=dist
            
            # 计算目标sv的覆盖度信息
            postprocessing.coverage(task.genotype_svs,task.lead_provider,config)

            #Determine genotypes for unmatched input SVs
            # 如果某个sv有覆盖度，则设置基因型为0，否则设置为none
            for svcall in task.genotype_svs:
                coverage_list=[svcall.coverage_start,svcall.coverage_center,svcall.coverage_end]
                coverage_list=[c for c in coverage_list if c!=None]
                if len(coverage_list)==0:
                    return # 不执行后续代码
                coverage=round(sum(coverage_list)/len(coverage_list)) # 平均覆盖度
                svcall.genotypes={}
                if coverage>0:
                    svcall.genotypes[0]=(0,0,0,coverage,0,None)
                else:
                    svcall.genotypes[0]=config.genotype_none

            result={}
            result["task_id"]=task.id
            result["processed_read_count"]=read_count
            result["svcalls"]=task.genotype_svs # 目标sv
            pipe.send(["return_genotype_vcf",result]) # 返回结果
            del task
            gc.collect()

        elif command=="combine":
            # combine模式
            task=arg
            result={}

            result["svcalls"],candidates_processed=task.combine(config)

            result["task_id"]=task.id
            result["processed_read_count"]=candidates_processed
            pipe.send(["return_combine",result])
            del task
            gc.collect()

        elif command=="finalize":
            # finalize模式，进程退出
            return
