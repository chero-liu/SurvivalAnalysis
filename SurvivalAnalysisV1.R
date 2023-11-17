library(survival)
library(survminer)
library(RTCGA.clinical)
library(ggplot2)
library(argparser)
source("/SurvivalAnalysis/plot.r")
# Command-line argument parsing
parser <- arg_parser(description = "Survival Analysis and Plotting Script")
parser <- add_argument(parser, "--intgene", help = "Comma-separated list of genes")
parser <- add_argument(parser, "--disease", help = "TCGA disease name")
parser <- add_argument(parser, "--input", help = "Input data file path")
parser <- add_argument(parser, "--outdir", help = "Output directory")
parser <- add_argument(parser, "--cutoff", help = "Cutoff value (e.g., median, best, 0.25)", default = "median")
parser <- add_argument(parser, "--prefix", help = "Output file prefix", default = "result")
parser <- add_argument(parser, "--timetype", help = "years,month,days", default = "days")
parser <- add_argument(parser, "--level", help = "intgene high or low", default = "None")
# Parse command-line arguments
args <- parse_args(parser)

# Get values from command-line arguments
intgene <- unlist(strsplit(args$intgene, ","))
disease <- args$disease
input <- args$input
outdir <- args$outdir
cutoff <- args$cutoff
prefix <- args$prefix
timetype <- args$timetype
level <- args$level

outdir = paste0(outdir,'/',disease)
dir.create(outdir)

# Get clinical data
var_name <- paste(disease, '.clinical', sep = '')
survivalTCGA(get(var_name)) -> survInfo

# Read gene expression data
exp <- t(read.csv(input, sep = '\t', check.names = FALSE))
bcr_patient_barcode <- row.names(exp)
exp <- as.data.frame(exp)
genes <- colnames(exp)
exp$bcr_patient_barcode <- bcr_patient_barcode

# Merge gene expression data and clinical data
merge_data <- merge(exp, survInfo, by = 'bcr_patient_barcode')

# Remove unnecessary data
rm(exp)
rm(survInfo)

# Check if intgene is valid
matching_genes <- intgene[intgene %in% genes]
non_matching_genes <- intgene[!intgene %in% genes]

if (length(matching_genes) == 0) {
  stop("None of the genes in intgene match any genes in the expression data. The program is terminated.")
}else if(length(non_matching_genes) > 0){
    cat(paste0(non_matching_genes, ' not in the expression data\n'))
}

if(level == 'None'){
    # Iterate over each matching gene for survival analysis and plotting
    for (gene in matching_genes){
        subdata = merge_data[c('times','patient.vital_status',gene)]
        # cutoff group
        if(cutoff=='median'){
        subdata$group = ifelse( subdata[,gene] <= median(subdata[,gene]), 'low', 'high')
        }else if (cutoff=='best'){
        cutpoint = surv_cutpoint(subdata, time=paste0(survtag,'.time'), event=survtag, variables=gene)$cutpoint[1,1]
        subdata$group = ifelse( subdata[,gene] <= cutpoint, 'low', 'high')
        }else{
        cut_value1=quantile(subdata[,gene], as.numeric(cutoff)) ;  cut_value2=quantile(subdata[,gene], (1-as.numeric(cutoff)) )
        print(c(cut_value1,cut_value2))
        subdata$group = ifelse(subdata[,gene] <= cut_value1,'low',ifelse(subdata[,gene] > cut_value2,'high','NA'))
        subdata = subdata[subdata$group %in% c('low','high'),]
        }
        
        #save table
        write.csv(subdata,paste0(outdir,'/',prefix,'_',cutoff,'_',gene,'.csv'))


        if(timetype == 'years'){
            subdata$times = subdata$times / 366
            xlab = 'Time (Years)'
        }else if(timetype == 'month'){
            subdata$times = subdata$times / 30
            xlab = 'Time (Month)'
        }else{
            xlab = 'Time (Day)'
        }
        #plot
        fit <- surv_fit( as.formula(Surv(times,patient.vital_status)~group), data=subdata)

        cph <- summary( coxph(Surv(times,patient.vital_status)~group,subdata) )
        line <- c( round(cph$coef[2],2), round(cph$conf.int[3:4],2),  round(cph$waldtest[3],2) )
        HR <- paste0(line[1],' (',line[2],'-',line[3],')')

        plot_ggsurvplot(fit,xlab = xlab,HR = HR,outdir = paste0(outdir,'/',prefix,'_',cutoff,'_',gene))
    }
}else{
    if(length(matching_genes) != 2){
        stop('Please sure matching_genes == 2')
    }
    subdata = merge_data[c('times','patient.vital_status',matching_genes)]
    geneA = matching_genes[1]
    geneB = matching_genes[2]

    res.cat = process_cutoff(data = subdata,
                geneA = geneA,
                geneB = geneB,
                outdir = outdir,
                cutoff = cutoff,
                prefix = prefix)

    table_name = paste0(outdir,'/',geneA,"_",geneB,"_",prefix,'_rescat.xls')
    write.table(res.cat, table_name, sep='\t', col.names=NA, row.names=T, quote=F)

    level1 <- unlist(strsplit(level,','))[1]
    level2 <- unlist(strsplit(level,','))[2]
    level3 <- unlist(strsplit(level,','))[3]
    level4 <- unlist(strsplit(level,','))[4]

    res.cat <- subset_group(res.cat,geneA,geneB,level1,level2,level3,level4)
    res.cat$group <- paste(res.cat[,geneA],geneA,"+",res.cat[,geneB],geneB)

    if(timetype == 'years'){
        res.cat$times = res.cat$times / 366
        xlab = 'Time (Years)'
    }else if(timetype == 'month'){
        res.cat$times = res.cat$times / 30
        xlab = 'Time (Month)'
    }else{
        xlab = 'Time (Day)'
    }

        
    fit <- surv_fit(Surv(times,patient.vital_status)~group,res.cat)

    cph <- summary( coxph(Surv(times,patient.vital_status)~group,res.cat) )
    line <- c( round(cph$coef[2],2), round(cph$conf.int[3:4],2),  round(cph$waldtest[3],2) )
    HR <- paste0(line[1],' (',line[2],'-',line[3],')')

    plot_ggsurvplot(fit,xlab = xlab,HR = HR,outdir = paste0(outdir,'/',geneA,'_',geneB,'_',cutoff,'_',prefix))
}

