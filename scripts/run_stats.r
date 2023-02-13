#!/usr/bin/Rscript

output_dir <- as.character(commandArgs(TRUE)[1])
mapping_rate <- as.character(commandArgs(TRUE)[2])

tbl <- read.table(paste0(output_dir,"/lengths_by_transcripts.txt"), sep = "\t")
colnames(tbl) <- c("tr_id","length","per_AT","per_CG","per_N","gene","tr_num")

num_gen <- length(unique(tbl$gene))
num_tr <- length(unique(tbl$tr_id))

tr_mayor500bp <- length(unique(tbl$tr_id[tbl$length>500]))
tr_mayor1000bp <- length(unique(tbl$tr_id[tbl$length>1000]))
gen_mayor500bp <- length(unique(tbl$gene[tbl$length>500]))
gen_mayor1000bp <- length(unique(tbl$gene[tbl$length>1000]))
per_AT <- round(mean(tbl$per_AT),2)
per_CG <- round(mean(tbl$per_CG),2)


#N50 y N90

N_X_L <- function(porcentaje, vector){
  total_esperado <- sum(vector)
  vector_ordenado <- vector[order(vector, decreasing = T)]
  NX_esperado <- total_esperado*(porcentaje/100)
  
  res_nx <- 0
  i_nx <- 1
  while(res_nx < NX_esperado) {
    x <- as.vector(vector_ordenado[i_nx])
    res_nx <- sum(res_nx+x)
    i_nx <- i_nx+1
    }
  
  NX <- i_nx
  LX <- vector_ordenado[i_nx]
  
  m <- c(NX,LX)
  names(m) <- c(paste0("N",porcentaje), paste0("L",porcentaje))
return(m)                
}

#Por Transcriptos:tomando isoformas para calcular el N50 y N90

N50_N_tr <- N_X_L(50,tbl$length)[1] # numero de trascriptos menor igual
N50_L_tr <- N_X_L(50,tbl$length)[2] #la longitud del transcripto

N90_N_tr <- N_X_L(90,tbl$length)[1] # numero de trascriptos menor igual
N90_L_tr <- N_X_L(90,tbl$length)[2] #la longitud del transcripto

#Por genes:tomando el transcripto de mayor longitud por cada gen
max_len_by_gene <- tapply(tbl$length, list(tbl$gene), max)

N50_N_gen <- N_X_L(50,max_len_by_gene)[1] # numero de trascriptos menor igual
N50_L_gen <- N_X_L(50,max_len_by_gene)[2] #la longitud del transcripto

N90_N_gen <- N_X_L(90,max_len_by_gene)[1] # numero de trascriptos menor igual
N90_L_gen <- N_X_L(90,max_len_by_gene)[2] #la longitud del transcripto


#Nx50 y Nx90

tabla_genes <- tbl[,c("tr_id","gene")]

files <- list.files(paste0(output_dir,"/salmon_quant"))

if(all(file.exists(files))){
  tpm <- tximport::tximport(files, type = "salmon",tx2gene = tabla_genes, countsFromAbundance = "scaledTPM")
  
}

tpm_sumados <- apply(tpm$abundance, 1,sum)
longitud_ponderada <- rowSums(tpm$length*tpm$abundance)/rowSums(tpm$abundance,)

df_expresion <- data.frame(Ex=NA, num_genes=NA,Expression_value=NA, ExN50=NA, ExpRel=NA)

for(i in seq(2,100,2)){
  res <- c(i,N_X_L(i,tpm_sumados))
  df_expresion <- rbind(df_expresion, res)
}


tpm_sumados <- tpm_sumados[order(tpm_sumados, decreasing = T)]

df_expresion <- data.frame()
for(i in 1:50){
  
  per_ex <- seq(2,100,2)[i]
  posicion <- which(as.numeric(tpm_sumados)==N_X_L(per_ex,tpm_sumados)[2])[1]-1
  tpm_x <- tpm_sumados[1:posicion]
  names_tpm_x <- names(tpm_x)
  subset_len <- longitud_ponderada[names(longitud_ponderada)%in%names_tpm_x]
  
  res <- data.frame(
  per_ex=per_ex,
  ExpRel=sum(tpm_x)/sum(tpm_sumados)*100,
  ExN50=N_X_L(50,subset_len)[2],
  per_genes=posicion/length(tpm_sumados)*100
  )
  
  df_expresion <- rbind(df_expresion,res)
}


pdf(file=paste0(output_dir,"/ExN50_plot.pdf"))

{
par(oma=c(5,5,5,5))
plot(df_expresion$per_ex, df_expresion$ExN50,type="b",pch=16, col="red", ylab="ExN50 (red)", xlab="Expression percentage")
abline(v=90)
par(new=TRUE)
plot(df_expresion$per_ex, df_expresion$per_genes,type="b",pch=16, col="blue", 
     ylab="", xlab="", axes=F)
axis(side = 4, at= pretty(range(df_expresion$per_genes)))
mtext("Gene percentage (blue)",side = 4, line = 3)
}

dev.off()


write.table(df_expresion,file = paste0(output_dir,"/ExN50_table.tsv"), sep="\t",quote = F,row.names = F)

num_gen <- length(unique(tbl$gene))
num_tr <- length(unique(tbl$tr_id))

tr_mayor500bp <- length(unique(tbl$tr_id[tbl$length>500]))
tr_mayor1000bp <- length(unique(tbl$tr_id[tbl$length>1000]))
gen_mayor500bp <- length(unique(tbl$gene[tbl$length>500]))
gen_mayor1000bp <- length(unique(tbl$gene[tbl$length>1000]))
per_AT <- round(mean(tbl$per_AT),2)
per_CG <- round(mean(tbl$per_CG),2)
mdian_tr_len <- median(tbl$length)
mdian_gene_len <- median(max_len_by_gene)
Ex90N50 <- df_expresion$ExN50[df_expresion$per_ex==90]
Gene_per_Ex90 <- df_expresion$per_genes[df_expresion$per_ex==90]


summary_table <- data.frame(
  Estatistics=c("Total 'genes'","Total transcripts","Genes >500bp (*)","Genes >1000bp (*)", 
                "Transcripts >500bp","Transcripts >1000bp","Percentage CG","Percentage AT",
    "N50 genes (*)", "N50L genes (*)","N90 genes (*)", "N90L genes (*)", "N50 transcripts", "N50L transcripts",
    "N90 transcripts", "N90L transcripts", "Median genes length (*)","Median transcripts length", 
    "Ex90N50 (**)", "Gene percentage Ex90", "Mapping rate"),
  
  Values=round(c(num_gen, num_tr, gen_mayor500bp, gen_mayor1000bp, tr_mayor500bp,tr_mayor1000bp, 
                 per_CG, per_AT, N50_N_gen, N50_L_gen, N90_N_gen, N90_L_gen, N50_N_tr, N50_L_tr, 
                 N90_N_tr, N90_L_tr, mdian_gene_len, mdian_tr_len,Ex90N50,Gene_per_Ex90, mapping_rate),2)

)

write.table(summary_table,file = paste0(output_dir,"/summary_table.tsv"), sep="\t",quote = F,row.names = F)
write.table(summary_table,file = paste0(output_dir,"/summary_table.csv"), sep=",",quote = F,row.names = F)

print(summary_table)
print('(*) the gene lengths is taken as the longest transcript for each gene')
print('(**) the gene length is taken as the expression-weighted mean of isoform lengths')

