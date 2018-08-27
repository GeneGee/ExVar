'''
Author: Zhang Zhao
Start Date:2017/8/30
Update: 2018/08/13
Version:2.0
Purpose: to give a variant a proper pathogenic/benign score, and a right class according to that score
Object: variant called (a line in the vcf file)
Class: Variant
Attribute: 
population_data; computation_data; function_data; segregation_data;
Method: 
pathogenic:pvs_score(pvs1); ps_score(ps1~4); pm_score(pm1~6); pp_score(pp1~5) 
benign:ba_score(ba1); bs_score(bs1~4); bp_score(bp1~7)
Specification: 
'''
import sys
#import gzip
import re
import os
import argparse

parser = argparse.ArgumentParser(prog='vep', description='VCE: Classify germline variants according to 2015 ACMG Guidelines',\
                                 epilog='Author: Zhang Zhao, Version: 2.0, email: zhangzhao1@genomics.cn')
parser.add_argument('mode')
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('log')

parser.add_argument('--lof',metavar='loss_of_function_intolerant',help='LOF intolerant genes',required=False)
parser.add_argument('--patho','-p', metavar='known_patho',help='pathogenic sites file',required=False)
parser.add_argument('--ps4',metavar='ps4_variants',help='obviously higher MAF in case than in control',required=False)
parser.add_argument('--pm1',metavar='pm1_sites',help='hotspot sites',required=False)
parser.add_argument('--iso',metavar='iso_splicing',help='splicing variants which result in in-frame RNA isomers',required=False)
parser.add_argument('--mis',metavar='mis_splicing',help='missense variants which have splicing effects',required=False)
parser.add_argument('--nifty',metavar='nifty_maf',help='population frequency in NIFTY',required=False)
parser.add_argument('--bp1',metavar='bp1_genes',help='genes whose truncations are the only known pathogenic mechanism',required=False)
parser.add_argument('--pp2',metavar='pp2_genes',help='genes in which missense variation is a common cause of disease ',required=False)
parser.add_argument('--exon',metavar='exon_info',help='number of exons and length of protein translated',required=False)
#parser.add_argument('--cancer',metavar='only_cancer',help='in this mode, we only classify genes related with cancer',required=False)
parser.add_argument('--morb',metavar='morbidity',help='user defined morbidity',default=0.001, type=float)
args = parser.parse_args()


def load_tumor_gene():
    gene_list = []
    with open(args.cancer) as fr_cancer:
        for line in fr_cancer.readlines():
            lineArr = line.strip().split()
            if lineArr[0] == 'Gene': continue
            gene_list.append(lineArr[0])
    return gene_list

def load_lof_data():
    lof_dict = {}
    with open(args.lof) as fr_lof:
        for line in fr_lof.readlines():
            lineArr = line.strip().split('\t')
            lof_dict[(lineArr[0],lineArr[1])] = [lineArr[2],lineArr[3]]
    return lof_dict
   
def load_bp1_data():
    bp1_list = []
    with open(args.bp1) as fr:
        for line in fr.readlines():
            bp1_list.append(line.strip())
            #print(bp1_list)
    return bp1_list
    
def load_ps4_data():
    ps4_list = []
    with open(args.ps4) as fr:
        for line in fr.readlines():
            lineArr = line.strip().split()
            ps4_list.append(','.join([lineArr[0],lineArr[1],lineArr[3],lineArr[4]]))          
    return ps4_list 
    
def load_pm1_data():
    pm1_dict = {}
    with open(args.pm1) as fr:
        for line in fr.readlines():
            lineArr = line.strip().split()
            key = lineArr[0]
            value = ','.join(lineArr[1:])
            if key not in pm1_dict:
                pm1_dict[key] = [value]
            else:
                pm1_dict[key].append(value)          
    return pm1_dict
    
def load_pp2_data():
    pp2_list = []
    with open(args.pp2) as fr:
        for line in fr.readlines():
            pp2_list.append(line.strip())
    return pp2_list    
   
def load_patho_data():
    same_aa_dict = {} 
    same_site_dict = {}
    with open(args.patho) as f_pat:
        for line in f_pat.readlines():
            lineArr = line.strip().split('\t')
            transcript = lineArr[0]
            pChange = lineArr[3]
            cChange = lineArr[2]
            pubID = lineArr[4]
            pPos = re.match(r'p\.([A-Z])(\d+)([A-Z])',pChange).group(2)
            pAlt = re.match(r'p\.([A-Z])(\d+)([A-Z])',pChange).group(3)
            cPos = re.match(r'c\.([A-Z])(\d+)([A-Z])',cChange).group(2)
            cAlt = re.match(r'c\.([A-Z])(\d+)([A-Z])',cChange).group(3)
            aa_key = (transcript, pChange)
            si_key = (transcript, pPos)
            same_aa_dict[aa_key] = [cChange,pubID]
            same_site_dict[si_key] = [pChange,pubID]
        #print('same_site_dict is',same_site_dict)
        #print('same_aa_dict is',same_aa_dict)
    return same_aa_dict,same_site_dict

def load_iso_data():
    iso_list = []
    with open(args.iso) as f_iso:
        for line in f_iso.readlines():
            if '#' in line:
                continue
            lineArr = line.strip().split('\t')
            iso_list.append([lineArr[0],lineArr[1]])
    return iso_list
  
def load_mis_data():
    mis_list = []
    with open(args.mis) as f_mis:
        for line in f_mis.readlines():
            if '#' in line:
                continue
            lineArr = line.strip().split('\t')
            mis_list.append([lineArr[0],lineArr[1]])
    return mis_list

def load_nifty_data():
    nifty_dict = {}
    with open(args.nifty) as f_nif:
        for line in f_nif.readlines():
            if '#' in line:
                continue
            lineArr = line.strip().split('\t')
            key = (lineArr[0].strip('chr'), lineArr[1])
            base = ['A','C','T','G']
            for i in range(4):
                nifty_dict[(key,base[i])] = lineArr[i+3]
    return nifty_dict

def load_exon_info():
    transcript_dict = {}
    with open(args.exon) as f_exon:
        for line in f_exon.readlines():
            if '#' in line:
                continue
            lineArr = line.strip().split('\t')
            key = lineArr[0]
            transcript_dict[key] = (lineArr[1],lineArr[2])
    return transcript_dict

def load_ano_data(file):#ann_list.append([pos_list,alter_list,var_list,freq_list,pred_list,spl_list,clin_list])
    pos = [Chr, Start, End] = [0,1,2]
    alter = [Ref, Alt] = [3,4]
    var =  [Func, Gene, GeneDetail, ExonicFunc, AAChange] = [5,6,7,8,9]
    pred = [SIFT_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, MutationTaster_pred, MutationAssessor_pred] = [32,35,38,44,47]
    spl = [dbscSNV_ADA, dbscSNV_RF, dpsi_max_tissue, dpsi_zscore ] = [96,97,98,99]
    freq = [ExAC_ALL, ExAC_EAS, ALL, EAS] = [18,21,12,15]
    clin = [CLINSIG, CLNDN, CLNREVSTAT] = [104,101,103]
    dom = [Interpro_domain] = [93]
    rmsk = [repetitive_region] = [105]
    ann_list = []
    pos_list = []
    alter_list = []
    var_list = []
    freq_list = []
    pred_list = []
    spl_list = []
    clin_list = []
    #dom_list = []
    with open(file,'rb') as f_ann:
        #print(f_ann)
        for line in f_ann.readlines():
            #print(line)
            lineArr = line.decode().strip().split('\t')
            if lineArr[0] == 'Chr':
                continue
            pos_list = [lineArr[pos[i]] for i in range(len(pos))]
            #print("pos_list is\n")
            #print(pos_list)
            alter_list = [lineArr[alter[i]] for i in range(len(alter))]
            #print("alter_list is\n")
            #print(alter_list)
            var_list = [lineArr[var[i]] for i in range(len(var))]
            #print("var_list is \n")
            #print(var_list)
            freq_list = [lineArr[freq[i]] for i in range(len(freq))]
            #print("freq_list is \n")
            #print(freq_list)
            pred_list = [lineArr[pred[i]] for i in range(len(pred))]
            #print("pred_list is \n")
            #print(pred_list)
            spl_list = [lineArr[spl[i]] for i in range(len(spl))]
            #print("spl_list is \n")
            #print(spl_list)
            clin_list = [lineArr[clin[i]] for i in range(len(clin))]
            #print("clin_list is\n")
            #print(clin_list)
            #dom_list = [lineArr[dom[i]] for i in range(len(dom))]
            rmsk_list = [lineArr[rmsk[i]] for i in range(len(rmsk))]
            ann_list.append([pos_list,alter_list,var_list,freq_list,pred_list,spl_list,clin_list,rmsk_list])
    return ann_list
    
class Variant(object):
    def __init__(self, ann_list, num, db_dict):
        
        self.ann = ann_list[num]
        self.nifty = db_dict['nifty']
        self.smaa = db_dict['same_aa']
        self.smsi = db_dict['same_site']
        self.iso = db_dict['iso']
        self.mis = db_dict['mis']
        #self.dom = dom_list
        self.exon = db_dict['exon']
        self.bp1List = db_dict['bp1']
        self.pp2List = db_dict['pp2']
        self.lof = db_dict['lof']
        self.pm1Dict = db_dict['pm1']
        self.ps4List = db_dict['ps4']
        self.log = []
        self.keyword = ''
        
        self.chr = self.ann[0][0]
        #print(self.chr)
        self.start = self.ann[0][1]
        #print("self.start is "+self.start)
        self.end = self.ann[0][2]
        #print("self.end is "+self.end)
        self.ref = self.ann[1][0]
        #print("self.ref is "+self.ref)
        self.alt = self.ann[1][1]
        #print("self.alt is "+self.alt)
        self.func = self.ann[2][0]
        #print("self.func is "+self.func)
        self.gene = self.ann[2][1]
        #print("self.gene is "+self.gene)
        self.genedetail = self.ann[2][2]
        #print("self.genedetail is "+self.genedetail)
        self.exonicfunc = self.ann[2][3]
        #print("self.exonicfunc is "+self.exonicfunc)
        self.aachange = self.ann[2][4]
        #print("self.aachange is " +self.aachange)
        self.maf = self.ann[3]
        #print(self.maf)
        #print("\n")
        self.pred = self.ann[4]
        #print(self.pred)
        #print("\n")
        self.spl = self.ann[5]
        #print(self.spl)
        #print("\n")
        self.clin = self.ann[6]
        #print(self.clin)
        #print("\n")
        #print(self.aachange)
        self.rmsk = self.ann[7]
        self.change_list = re.split(r'[;,]', self.aachange)
        self.geneDetail_list = re.split(r'[;,]', self.genedetail)
        self.stop_pattern = r'p\.([A-Z])(\d+)([A-Z]fs\*\d+)?(\*|X)?'
        self.splicing_pattern = r'c\.(\d+)(\+|\-)(1|2)'
        self.aa_pattern = r'p\.([A-Z])(\d+)([A-Z])'
        self.cod_pattern = r'c\.([A-Z])(\d+)([A-Z])'
        self.patho = re.compile(r'pathogenic', re.I)
        self.benign = re.compile(r'benign', re.I)
        if self.func == 'exonic':
            #print(self.change_list)
            for item in self.change_list:
                if item == 'UNKNOWN':
                    self.log.append('该变异是发生在%s基因上的未知类型(UNKNOWN)变异'%(self.gene))
                    continue
                item_list = item.split(':')
                self.keyword = item_list[-1]
                self.transcript = item_list[1]
                #print(self.transcript)
                #print(item)
                ac = item
                exon_number = item_list[2].strip('exon')
                if self.exonicfunc == 'nonsynonymous SNV':
                    mutype = '错义变异(nonsynonymous SNV)'
                elif self.exonicfunc == 'synonymous SNV':
                    mutype = '同义变异(synonymous SNV)'
                elif self.exonicfunc == 'frameshift deletion':
                    mutype = '移码缺失(frameshift deletion)'
                elif self.exonicfunc == 'frameshift insertion':
                    mutype = '移码插入(frameshift insertion)'
                elif self.exonicfunc == 'frameshift substitution':
                    mutype = '移码替换(frameshift substitution)'
                elif self.exonicfunc == 'nonframeshift substitution':
                    mutype = '非移码替换(nonframeshift substitution)'
                elif self.exonicfunc == 'nonframeshift deletion':
                    mutype = '非移码缺失(nonframeshift deletion)'
                elif self.exonicfunc == 'nonframeshift insertion':
                    mutype = '非移码插入(nonframeshift insertion)'
                elif self.exonicfunc == 'stopgain':
                    mutype = '无义变异(stopgain)'
                elif self.exonicfunc == 'stoploss':
                    mutype = '终止密码子缺失变异(stoploss)'
                elif self.exonicfunc == 'unknown':
                    mutype = '未知类型变异(unknown)'
                #else:
                    #print(ann_list[num])
                #print(mutype)
                if self.transcript in self.exon:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的%s，该转录本共有%s个外显子，对应蛋白长度为%s个氨基酸'%(ac,self.transcript,exon_number,mutype,self.exon[self.transcript][0],self.exon[self.transcript][1]))
                else:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的%s'%(ac,self.transcript,exon_number,mutype))
                               
        elif self.func == 'splicing' and self.geneDetail_list != ['.']:
            for item in self.geneDetail_list:
                item_list = item.split(':')
                self.transcript = item_list[0]
                #print(item)
                exon_number = item_list[1].strip('exon')
                self.keyword = item_list[2]
                ac = item
                if (self.gene,self.transcript) in self.exon:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的剪接变异，该基因共有%s个外显子，对应蛋白长度为%s个氨基酸'%(ac,self.transcript,exon_number,self.exon[self.transcript][0],self.exon[self.transcript][1]))
                else:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的剪接变异'%(ac,self.transcript,exon_number))
        elif self.func == 'intronic':
            self.log.append('该变异位于%s基因的内含子区域'%self.gene)
        elif self.func == 'exonic;splicing':
            for item in self.change_list:
                if item == 'UNKNOWN':
                    self.log.append('该变异是发生在%s基因上的未知类型(UNKNOWN)变异'%(self.gene))
                    continue
                item_list = item.split(':')
                self.keyword = item_list[-1]
                #print(item_list)
                self.transcript = item_list[1]
                
                #print(self.transcript)
                #print(item)
                ac = item
                exon_number = item_list[2].strip('exon')
                if self.exonicfunc == 'nonsynonymous SNV':
                    mutype = '错义变异(nonsynonymous SNV)'
                elif self.exonicfunc == 'synonymous SNV':
                    mutype = '同义变异(synonymous SNV)'
                elif self.exonicfunc == 'frameshift deletion':
                    mutype = '移码缺失(frameshift deletion)'
                elif self.exonicfunc == 'frameshift insertion':
                    mutype = '移码插入(frameshift insertion)'
                elif self.exonicfunc == 'frameshift substitution':
                    mutype = '移码替换(frameshift substitution)'
                elif self.exonicfunc == 'nonframeshift deletion':
                    mutype = '非移码缺失(nonframeshift deletion)'
                elif self.exonicfunc == 'nonframeshift insertion':
                    mutype = '非移码插入(nonframeshift insertion)'
                elif self.exonicfunc == 'stopgain':
                    mutype = '无义变异(stopgain)'
                elif self.exonicfunc == 'stoploss':
                    mutype = '终止密码子缺失变异(stoploss)'
                #print(mutype)
                if (self.gene,self.transcript) in self.exon:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的%s，该转录本共有%s个外显子，对应蛋白长度为%s个氨基酸'%(ac,self.transcript,exon_number,mutype,self.exon[self.transcript][0],self.exon[self.transcript][1]))
                else:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的%s'%(ac,self.transcript,exon_number,mutype))
            
            for item in self.geneDetail_list:
                item_list = item.split(':')
                #print(item_list)
                if item_list == ['.']: continue
                self.transcript = item_list[0]
                #print(item)
                exon_number = item_list[1].strip('exon')
                self.keyword = item_list[2]
                ac = item
                if (self.gene,self.transcript) in self.exon:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的剪接变异，该基因共有%s个外显子，对应蛋白长度为%s个氨基酸'%(ac,self.transcript,exon_number,self.exon[self.transcript][0],self.exon[self.transcript][1]))
                else:
                    self.log.append('%s是发生在%s转录本第%s号外显子上的剪接变异'%(ac,self.transcript,exon_number))
            
    def out_scope(self):
        is_out = 0
        if self.func not in ['exonic','splicing','exonic;splicing','intronic']:
            is_out = 1
        return is_out
        
    def pvs1(self): # wholegene deletion, stopgain, frameshift, or splicing mutation in 2bp 
        is_pvs1 = 0
        is_lof = 0
        for key in self.lof.keys():
            if self.gene in key:
                is_lof = 1
                lof_score = self.lof[key][0]
                lof_pos = int(self.lof[key][1])
        if is_lof == 0:
            return is_pvs1
        else:
            self.log.append('其LoFtool分值为%s'%lof_score)
        if self.func == 'exonic':
            if self.exonicfunc == 'stopgain' or ('frameshift' in self.exonicfunc and 'nonframeshift' not in self.exonicfunc):
                #print(self.exonicfunc)
                pos = 0
                item = self.change_list[-1]
                item_list = item.split(':')
                #print(item_list)
                if self.aachange == '.':
                    self.log.append('未成功预测氨基酸发生的改变(PVS1)') #no amino acid change predicted
                    return is_pvs1
                elif item_list[2] == 'wholegene':
                    is_pvs1 = 1
                    self.log.append('起始密码子丢失导致了%s基因完全的丢失，构成PVS1证据'%self.gene) #'PVS1:全基因缺失':
                    is_pvs1 = 1
                elif len(item_list) == 4:
                    transcript = item_list[1]
                    if (self.gene,transcript) in self.lof:
                        print(self.exonicfunc)
                        print(item_list)                    
                        is_pvs1 = 1
                        self.log.append('移码替换，构成PVS1证据') #'PVS1:全基因缺失':
                        return is_pvs1

                elif re.match(self.stop_pattern, item_list[4]):
                    pos = int(re.match(self.stop_pattern, item_list[4]).group(2))
                    transcript = item_list[1]
                    if (self.gene,transcript) not in self.lof: return is_pvs1
                    lof_pos = int(self.lof[(self.gene,transcript)][1])
                    if lof_pos > pos:
                        is_pvs1 = 1
                        self.log.append('截短变异在第%d位氨基酸处，构成PVS1证据'%pos)
                    return is_pvs1
                    
        elif self.func == 'splicing':
            nam = ''
            judge = 1
            if self.genedetail == '.':
                return is_pvs1
            for item in self.geneDetail_list:
                item_list = item.split(':')
                if (item_list[0] == 'NM_007294') or (item_list[0] == 'NM_000059'):
                    if re.match(self.splicing_pattern, item_list[2]):
                        nam = re.match(self.splicing_pattern, item_list[2]).group(0)
                        for iso in self.iso:
                            if iso[0] == item_list and iso[1] in nam: 
                                judge = 0
                                self.log.append('该变异是经典2bp以内的剪接变异，但会导致框内RNA异构体(不构成PVS1证据)')
            if judge == 1:                
                is_pvs1 = 1
                self.log.append('该变异是经典2bp以内的剪接变异(PVS1)')
                return is_pvs1
                #PVS1: canonical splicing variant which does not result in in-frame RNA isomer
        return is_pvs1
    
    def ps1(self): # same amino acid change as known pathogenic variants, with different nucleic acid alteration 
        is_ps1 = 0
        sep = '\t'
        judge = 0
        if self.func == 'exonic':           
            if self.exonicfunc == 'nonsynonymous SNV':
                for aa in self.change_list:
                    item_list = aa.split(':')
                    key = (item_list[1], item_list[4])
                    value = item_list[3]
                    if key in self.smaa and value != self.smaa[key][0]:
                        pub_link = self.smaa[key][1]
                        judge = 1
                        break
        if judge == 1:
            is_ps1 = 1
            self.log.append('该变异与已知的具有致病性的错义变异在同一位点具有相同的氨基酸改变，但碱基改变不同(PS1),参考文献%s'%pub_link)
            #'PS1: the same amino acid change as a known pathogenic missense variant, with different nucleic acid change
        return is_ps1
        
    def ps2(self): # de novo variant 
        is_ps2 = 0
        return is_ps2
        
    def ps3(self):# destructive effects verified by in vitro/vivo function tests
        is_ps3 = 0
        return is_ps3
    
    def ps4(self):# obviously higher frequency in case than in control 
        is_ps4 = 0
        key = ','.join([self.chr, self.start, self.ref, self.alt])
        if key in self.ps4List:
            is_ps4 = 1
        return is_ps4
        
    def pm1(self): # located in the hotspot area or specific function domain 
        is_pm1 = 0
        if self.gene in self.pm1Dict:
            if [self.chr,self.start] in self.pm1Dict[self.gene]:
                is_pm1 = 1
        return is_pm1
        
    def pm2(self): # no record OR lower than 0.0003 in 1000g2015aug(ALL,EAS), exac03(ALL,_EAS)
        is_pm2 = 0
        freq_list = self.maf
        judge1, judge2 = 1, 1
        freq_nifty = 0.0
        for fq in freq_list:
            if fq != '.' and float(fq) > 0.0003:
                judge1 = 0
                '''
        if 'SNV' in self.exonicfunc:
            key = ((self.chr, self.start),self.alt) 
            if key in self.nifty:
                if float(self.nifty[key]) > 0.0003:
                    judge2 = 0
                else:
                    freq_nifty = self.nifty[key]
                    '''
        if judge1 == 1:# and judge2 == 1:
            is_pm2 = 1
            self.log.append('该变异为稀有变异(PM2)')
            #'PM2:在ExAC、千人基因组和NIFTY数据库中的频率低于万分之三'
        return is_pm2
    
    def pm3(self): # detected in trans in recessive genetic disorder
        is_pm3 = 0
        return is_pm3 
        
#An in-frame deletion/insertion is a mutation in the coding region of a gene, resulting in the loss/gain of DNA sequence, which does not alter the normal triplet reading frame between the mutation site and the carboxyterminus of the polypeptide     
    def pm4(self): # in-frame deletion/insertion OR stop codon loss which resulted the length change of translated protein
        is_pm4 = 0
        if self.func == 'exonic' and self.rmsk == '.':
            if self.exonicfunc == 'stoploss':
                is_pm4 =1
                self.log.append('终止密码子缺失导致翻译出来的蛋白质长度发生了改变(PM4)')
                #PM4:stop codon loss which resulted the length change of translated protein
            elif 'nonframeshift' in self.exonicfunc:
                is_pm4 =1
                if 'insertion' in self.exonicfunc:
                    self.log.append('非重复区域的框内插入导致翻译出来的蛋白质长度发生了改变(PM4)')
                #PM4:in-frame deletion/insertion which resulted the length change of translated protein
                elif 'deletion' in self.exonicfunc:
                    self.log.append('非重复区域的框内缺失导致翻译出来的蛋白质长度发生了改变(PM4)')
        return is_pm4
        
    def pm5(self):# same site as the known pathogenic variant with different nucleic acid change
        is_pm5 = 0
        judge = 0
        if self.func == 'exonic':
            if self.exonicfunc == 'nonsynonymous SNV':
                for aa in self.change_list:
                    item_list = aa.split(':')
                    pos = re.match(r'p\.([A-Z])(\d+)([A-Z])',item_list[4]).group(2)
                    key = (item_list[1],pos)
                    value = item_list[4]
                    if key in self.smsi and value != self.smsi[key][0]:
                        pub_link = self.smsi[key][1]
                        judge = 1
                        break
        if judge == 1:
            is_pm5 = 1
            #print('PM5:a novel missense amino acid change ocurring at the same position as another pathogenic missense change')
            self.log.append('新的错义变异导致氨基酸变化，此变异之前未曾报道，但是在同一位点，导致另外一种氨基酸的变异已经确认是致病性的(PM5)，参考文献%s'%pub_link)
            #PM5:a novel missense amino acid change ocurring at the same position as another pathogenic missense change
        return is_pm5
      
    def pm6(self):# Assumed de novo, but without confirmation of paternity and maternity 
        is_pm6 = 0
        return is_pm6
      
    def pp1(self):# cosegregation with disease in multiple affected family members
        is_pp1 = 0
        return is_pp1
       
    def pp2(self):# missense variant in a gene that has a low rate of benign missense variation
        is_pp2 = 0
        if self.gene not in self.pp2List:
            return is_pp2
        if self.func == 'exonic':
            if self.exonicfunc == 'nonsynonymous SNV':
                is_pp2 = 1
                self.log.append('已知%s基因上大部分的致病变异为错义变异，且良性错义变异概率极低，而该变异恰为错义变异(PP2)'%self.gene)
                #missense variant while most pathogenic variants are truncating variants
        return is_pp2 

#SIFT:  D: Deleterious (sift<=0.05); T: Tolerated (sift>0.05)
#PolyPhen 2 HDIV: D: Probably damaging (>=0.957), P: Possibly damaging (0.453<=pp2_hdiv<=0.956), B: Benign (<=0.452)
#PolyPhen 2 HVar: D: Probably damaging (>=0.909), P: Possibly damaging (0.447<=pp2_hdvar<=0.909), B: Benign (<=0.446)
#MutationTaster: A: Disease_causing_automatic, D: Disease_causing, N: Polymorphism, P: Polymorphism_automatic
#MutationAssessor: H: High; M:Medium; L:Low, N:Neutral, H/M: functional, L/N: non-functional         
    def pp3(self):# multiple lines of computational evidence support a deleterious effect on the gene or gene product
        is_pp3 = 0
        SIFT = self.pred[0]
        PolyPhen_hdiv = self.pred[1]
        PolyPhen_hvar = self.pred[2]
        Taster = self.pred[3]
        Assessor = self.pred[4]
        if SIFT == 'D' and PolyPhen_hdiv == 'D' and PolyPhen_hvar == 'D':
            if Taster == 'A' or Taster == 'D': 
                if Assessor == 'H' or Assessor == 'M':
                    is_pp3 = 1
                    self.log.append('多种软件（SIFT/PolyPhen/MutationTaster/MutationAssessor）预测该变异会对基因或基因产物造成有害影响(PP3)')
                    #multiple lines of computational evidence support a deleterious effect on the gene or gene product, \
                    #SIFT: %s, PolyPhen_hdiv: %s, PolyPhen_hvar: %s, MutationTaster: %s, MutationAssessor: %s' %(self.pred[0],self.pred[1],self.pred[2],self.pred[3],self.pred[4]))
        return is_pp3
        
    def pp4(self):# patient's phenotype or family history is highly specific for a disease wiht a single genetic etiology
        is_pp4 = 0
        return is_pp4
        
    def pp5(self):# reputable source recently reports variant as pathogenic  clinvar/BIC
        is_pp5 = 0
        jud = self.clin[0]
        jud_list = re.split(r'[\|,]',jud)
        for it in jud_list:
            if self.patho.search(it):
                is_pp5 = 1
                self.log.append('有报道(ClinVar)认为该变异为具有致病性(PP5)，相关表型为%s，具体报道情况：%s'%(self.clin[1],self.clin[2]))
                #PP5:reputable source recently reports variant as pathogenic
                break
        return is_pp5
                
    def ba1(self):# Allele frequency is >1% 
        is_ba1 = 0
        freq_list = self.maf
        judge1, judge2 = 0, 1
        for fq in freq_list:
            if fq == '.':
                continue
            if float(fq) > 0.01:
                judge1 = 1
                '''
        if 'SNV' in self.exonicfunc:
            key = ((self.chr, self.start),self.alt) 
            if key in self.nifty:
                if float(self.nifty[key]) <= 0.01:
                    judge2 = 0
                    '''
        if judge1 == 1:# and judge2 == 1:
            is_ba1 = 1
            self.log.append('在ExAC、千人基因组数据库中的频率均高于百分之一(BA1)')
            #frequencies in ExAC,1kG,and NIFTY database higher than 0.01
        return is_ba1
    
    def bs1(self):# Allele frequency is greater than expected for disorder
        is_bs1 = 0
        freq_list = self.maf
        judge1 = 0
        #print(args.morb)
        for fq in freq_list:
            if fq == '.':
                continue
            if float(fq) >= args.morb: #frequency of disorder
                judge1 = 1
                    
        if judge1 == 1: 
            is_bs1 = 1
            self.log.append('在ExAC、千人基因组数据库中的频率均高于发病率(BS1)')
        return is_bs1
       
    def bs2(self):
        is_bs2 = 0
        return is_bs2
        
    def bs3(self):
        is_bs3 = 0
        return is_bs3
        
    def bs4(self):
        is_bs4 = 0
        return is_bs4
  
    def bp1(self):# missense variants
        is_bp1 = 0
        if self.gene not in self.bp1List:
            return is_bp1
        if self.func == 'exonic':
            if self.exonicfunc == 'nonsynonymous SNV':
                is_bp1 = 1
                self.log.append('已知%s基因上大部分的致病变异为截短变异，而该变异为错义变异(BP1)'%self.gene)
                #missense variant while most pathogenic variants are truncating variants
        return is_bp1
    
    def bp2(self):
        is_bp2 = 0
        return is_bp2 
    
    def bp3(self):# in-frame deletions/insertions in a repetitive region without a known function
        is_bp3 = 0
        if self.func == 'exonic' and self.rmsk != '.':
            if 'nonframeshift' in self.exonicfunc:
                is_bp3 =1
                if 'insertion' in self.exonicfunc:
                    self.log.append('重复区域的框内插入导致翻译出来的蛋白质长度发生了改变(BP3)')
                #PM4:in-frame deletion/insertion which resulted the length change of translated protein
                elif 'deletion' in self.exonicfunc:
                    self.log.append('重复区域的框内缺失导致翻译出来的蛋白质长度发生了改变(BP3)')
        return is_bp3
        
    def bp4(self):# multiple lines of computational evidence suggest no impact on gene or gene product
        is_bp4 = 0
        SIFT = self.pred[0]
        PolyPhen_hdiv = self.pred[1]
        PolyPhen_hvar = self.pred[2]
        Taster = self.pred[3]
        Assessor = self.pred[4]
        if SIFT == 'T':
            if PolyPhen_hdiv == 'B':
                if PolyPhen_hvar == 'B':
                    if Taster == 'N' or Taster == 'P': 
                        if Assessor == 'L' or Assessor == 'N':
                            is_bp4 = 1
                            self.log.append('多种软件（SIFT/PolyPhen/MutationTaster/MutationAssessor）均预测该变异对基因或基因产物无有害影响(BP4)')
                            #multiple lines of computational evidence suggest no impact on gene or gene product
        return is_bp4
        
    def bp5(self):
        is_bp5 = 0
        return is_bp5 
    
    def bp6(self):# reputable source recently reports variant as benign
        is_bp6 = 0
        jud = self.clin[0]
        jud_list = re.split(r'[\|,]',jud)
        for it in jud_list:
            if self.benign.search(it):
                is_bp6 = 1
                self.log.append('有报道(ClinVar)认为该变异为良性或可能良性(BP6)')
                #reputable source recently reports variant as benign
                break
        return is_bp6
       
    def bp7(self):# a synonymous variant for which splicing prediction algorithms predict no impact to the splice consensus sequence 
        is_bp7 = 0
        if self.exonicfunc == 'synonymous SNV':#nor the creation of a new splice site AND the nucleotide is not highly conserved
            if self.spl[0] == '.' or self.spl[1] == '.':
                is_bp7 = 1
                self.log.append('同义变异，且预测不影响剪接(BP7)')
                return is_bp7
            if float(self.spl[0]) <= 0.6 and float(self.spl[1]) <= 0.6:
                is_bp7 = 1
                self.log.append('同义变异，且预测不影响剪接(BP7)')
                # a synonymous variant for which splicing prediction algorithms predict no impact to the splice consensus sequence
        return is_bp7
        
    def log_wr(self):
        return self.keyword, self.log

def pvs_score(file, db_dict):
    title = ['sample','gene','var_type','chr','start','end','ref','alt', 'dbSNPID','amino acid change','cDNA change',\
            'class','summary','explanation','tips',\
            'pvs1','ps1','ps2','ps3','ps4','pm1','pm2','pm3','pm4','pm5','pm6','pp1','pp2','pp3','pp4',\
            'pp5', 'ba1','bs1','bs2','bs3','bs4','bp1','bp2','bp3','bp4','bp5','bp6','bp7']
    ann_list = load_ano_data(file)
    n = len(ann_list)
    var_class_list = []
    for i in range(n):
        acmg_class = ''
        line = []
        tip = []
        
        p = 0
        lp = 0
        b = 0
        lb = 0
        us = 0
        
        pvs = [0]
        ps = [0,0,0,0]
        pm = [0,0,0,0,0,0]
        pp = [0,0,0,0,0]
        ba = [0]
        bs = [0,0,0,0]
        bp = [0,0,0,0,0,0,0]
        
        var_in = Variant(ann_list,i,db_dict)
        if var_in.out_scope():
            var_class_list.append(['.']*(len(title)-10))
            continue
            
        pvs[0] = var_in.pvs1()
        
        ps[0] = var_in.ps1()
        ps[1] = var_in.ps2()
        ps[2] = var_in.ps3()
        ps[3] = var_in.ps4()
        
        pm[0] = var_in.pm1()
        pm[1] = var_in.pm2()
        pm[2] = var_in.pm3()
        pm[3] = var_in.pm4()
        pm[4] = var_in.pm5()
        pm[5] = var_in.pm6()
        
        pp[0] = var_in.pp1()
        pp[1] = var_in.pp2()
        pp[2] = var_in.pp3()
        pp[3] = var_in.pp4()
        pp[4] = var_in.pp5()

        
        ba[0] = var_in.ba1()
        
        bs[0] = var_in.bs1()
        bs[1] = var_in.bs2()
        bs[2] = var_in.bs3()
        bs[3] = var_in.bs4()
        
        bp[0] = var_in.bp1()
        bp[1] = var_in.bp2()
        bp[2] = var_in.bp3()
        bp[3] = var_in.bp4()
        bp[4] = var_in.bp5()
        bp[5] = var_in.bp6()
        bp[6] = var_in.bp7()
        con = []
        for item in [pvs,ps,pm,pp,ba,bs,bp]:
            con.extend(item)
        for item in con:
            line.append(str(item))
            
        keywr, log_info_list = var_in.log_wr()
        if log_info_list:
            log_info = ';'.join(log_info_list)
        
        if pp[4] == 1 and bp[5] == 1:
            tip.append("在clinvar中既有致病又有良性（或可能良性）的报道，请根据具体情况作出相应判读https://www.ncbi.nlm.nih.gov/clinvar/")
        
         
        if sum(pvs) == 1:
            if sum(ps) >= 1 or sum(pm) >= 2 or (sum(pm)>=1 and sum(pp)>=1):
                p = 1
            elif sum(pp) >= 2:
                p = 1
        elif sum(ps) >= 2:
            p = 1
        elif sum(ps) == 1:
            if sum(pm) >= 3:
                p = 1
            elif sum(pm) >= 2 and sum(pp) >= 2:
                p = 1
            elif sum(pm) == 1 and sum(pp) >=4:
                p = 1
        if sum(pvs) == 1 and sum(pm) == 1:
            lp = 1
        elif sum(ps) == 1 and sum(pm) in [1,2]:
            lp = 1
        elif sum(ps) == 1 and sum(pp) >= 2:
            lp = 1
        elif sum(pm) >= 3:
            lp = 1
        elif sum(pm) == 2 and sum(pp) >= 2:
            lp = 1
        elif sum(pm) == 1 and sum(pp) >= 4:
            lp = 1
        if sum(ba) == 1 or sum(bs) >= 2:
            b = 1
        if sum(bs) == 1 and sum(bp) == 1:
            lb = 1
        elif sum(bp) >= 2:
            lb = 1
        if p == 1 and b == 0:
            acmg_class = 'pathogenic'
        elif lp == 1 and p == 0:
            acmg_class = 'likely pathogenic'
        elif b == 1 and p == 0:
            acmg_class = 'benign'
        elif lb == 1 and b == 0:
            acmg_class = 'likely benign'
        else:
            acmg_class = 'uncertain significance'
        if acmg_class == 'uncertain significance':
            #print(keywr)
            tip.append('本系统不能给出明确判读，请参考谷歌学术https://scholar.google.com.hk/scholar?hl=zh-CN&as_sdt=0%2C5&'+keywr+'&btnG='+'，或PubMed数据库'+'https://www.ncbi.nlm.nih.gov/pubmed/?term='+keywr+'%3EC'+'作出相应判读，若体内/体外功能实验已明确会导致基因功能受损，则可加入PS3证据')
        
#title = ['sample','gene','var_type','chr','start','end','ref','alt', 'amino acid change','cDNA change','class','summary','explanation','tips','pvs1','ps1','ps2','ps3','ps4','pm1','pm2','pm3','pm4','pm5','pm6','pp1','pp2','pp3','pp4',\
           # 'pp5', 'pp6', 'ba1','bs1','bs2','bs3','bs4','bp1','bp2','bp3','bp4','bp5','bp6','bp7']        
        summary = []
        #print(len(line))
        #print(line)
        #print(title)
        #print(line)
        for j in range(len(line)):
            ite = line[j]
            if ite == '1':
                summary.append(title[j+15])
        smay = '；'.join(summary)
        
        line.insert(0,acmg_class)
        line.insert(1,smay)
        line.insert(2,log_info)
        line.insert(3,'，'.join(tip))
        #print(acmg_class)
        var_class_list.append(line)
    return var_class_list

def main():
    if len(sys.argv) == 1:
        os.system("python %s -h" %sys.argv[0])
    sample_name = os.path.basename(args.input).split('.')[0]
    print("The sample name is %s" %sample_name)
    if not os.path.exists(sample_name):
        os.system('mkdir %s' %sample_name)
        
    fw_result = open(sample_name+'/'+args.output,'w', encoding = 'gb18030')
    fw_log = open(sample_name+'/'+args.log,'w')

    if args.mode == 'ann':
        input_file = args.input
    elif args.mode == 'vcf':
        os.system('perl /hwfssz1/ST_CANCER/POL/USER/zhangzhao1/variant_classifier/annovar/acmg_class/table_annovar.pl %s -vcfinput /hwfssz1/ST_CANCER/POL/USER/zhangzhao1/variant_classifier/annovar/humandb \
                -buildver hg19 -out %s/%s -remove\
                -protocol refGene,cytoBand,popfreq_all_20150413,dbnsfp33a,dbscsnv11,spidex,clinvar_20180603,rmsk,ensGene,avsnp150,icgc21 -operation g,r,f,f,f,f,f,r,g,f,f -nastring . -polish'\
                %(args.input,sample_name,sample_name))
        input_file = '%s/%s.hg19_multianno.txt' %(sample_name,sample_name)
    elif args.mode == 'avinput':
        os.system('perl /hwfssz1/ST_CANCER/POL/USER/zhangzhao1/variant_classifier/annovar/acmg_class/table_annovar.pl %s /hwfssz1/ST_CANCER/POL/USER/zhangzhao1/variant_classifier/annovar/humandb \
                -buildver hg19 -out %s/%s -remove\
                -protocol refGene,cytoBand,popfreq_all_20150413,dbnsfp33a,dbscsnv11,spidex,clinvar_20180603,rmsk,ensGene,avsnp150,icgc21 -operation g,r,f,f,f,f,f,r,g,f,f -nastring . -polish'\
                %(args.input,sample_name,sample_name))
        input_file = '%s/%s.hg19_multianno.txt' %(sample_name,sample_name)
    #ann_list = load_ano_data(input_file)
    db_dict = {}
    same_aa_dict,same_site_dict = load_patho_data()
    db_dict['same_aa'] = same_aa_dict
    db_dict['same_site'] = same_site_dict
    db_dict['iso'] = iso_list = load_iso_data()
    db_dict['mis'] = mis_list = load_mis_data()
    db_dict['nifty'] = nifty_dict = load_nifty_data()
    db_dict['lof'] = lof_dict = load_lof_data()
    db_dict['bp1'] = bp1_list = load_bp1_data()
    db_dict['pp2'] = pp2_list = load_pp2_data()
    db_dict['pm1'] = pm1_dict = load_pm1_data()
    db_dict['ps4'] = ps4_list = load_ps4_data()
    db_dict['exon'] = exon_dict = load_exon_info() 

    #p,lp,b,lb,us    
    var_class = pvs_score(input_file, db_dict)
    #print("there are %d var" %len(var_class))
    #print(var_class)
    with open(input_file,'rb') as fr_input:
        fr_ann = fr_input.readlines()

    #print("there are %d fr_ann" %num_var)
    title = ['sample','gene','var_type','chr','start','end','ref','alt', 'dbSNPID','amino acid change','cDNA change',\
            'class','summary','explanation','tips',\
            'pvs1','ps1','ps2','ps3','ps4','pm1','pm2','pm3','pm4','pm5','pm6','pp1','pp2','pp3','pp4',\
            'pp5', 'ba1','bs1','bs2','bs3','bs4','bp1','bp2','bp3','bp4','bp5','bp6','bp7']
    fw_result.write('\t'.join(title)+'\t\n')
    fw_log.write('\t'.join(fr_ann[0].decode().strip().split('\t'))+'\t\n')
    for i in range(1,len(var_class)+1):
        #print(var_class[i].encode("utf-8"))
        annArr = fr_ann[i].decode().strip().split('\t')
        aacha = re.split(r'[,;]',annArr[9])
        pchange = '.'
        cchange = '.'
        if '.' not in aacha and 'UNKNOWN' not in aacha :
            ah = aacha[-1]
            pchange = ah.split(':')[-1]
            cchange = ah.split(':')[-2]
        #print(cchange)
        
        var_type = ''
        if annArr[8] != '.':
            var_type = annArr[5]+','+annArr[8]
        else:
            var_type = annArr[5]
        #print(sample_name)
        #print(sample_name+'\t'+annArr[6]+'\t'+var_type+'\t'+'\t'.join([annArr[i] for i in [0,1,2,3,4]])+'\t'+pchange+'\t'+cchange+'\t'+'\t'.join(var_class[i-1])+'\t\n')
        #fw_log.write(fr_ann[i].decode())
        fw_result.write(sample_name+'\t'+annArr[6]+'\t'+var_type+'\t'+'\t'.join([annArr[i] for i in [0,1,2,3,4]])+'\t'+annArr[111]+'\t'+pchange+'\t'+cchange+'\t'+'\t'.join(var_class[i-1])+'\t\n')
    
if __name__ == "__main__":
    main()
