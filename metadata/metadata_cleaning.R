rm(list = ls())

library(pacman)
p_load(ShortRead, Biostrings, tidyverse, XML)

# 读取fastq前先在Linux解压fastq.gz文件

# 提取RJEB22874 - Oral 16s的barcodes，共86个patients
if (F) {
  setwd("~/Wargo/PRJEB22874/ERR2155099")
  
  buccal <- readFastq("buccal.fastq")
  
  buccal@id
  
  buccal_barcode_list <- buccal@id %>% 
    as.character() %>% 
    str_replace("(Wargo\\.[0-9]+)_.*new_bc=([A-Z]+)\\s.*", 
              "\\1\\ \\2") %>% 
    unique() %>% 
    as.data.frame()
  
  write.table(buccal_barcode_list, file = "PRJEB22874_barcode_list.txt", 
            quote = F, row.names = F, col.names = F, sep = "\n")
}

# 提取PRJEB22894 - Fecal 16s的barcodes，共43个patients
if (F) {
  setwd("~/Wargo/PRJEB22894/ERR2162225")
  
  fecal_16S <- readFastq("fecal_16S.fastq")
  
  fecal_16S@id
  
  fecal_16S_barcode_list <- fecal_16S@id %>% 
    as.character() %>% 
    str_replace("(Wargo\\.[0-9]+)_.*new_bc=([A-Z]+)\\s.*", 
                "\\1\\ \\2") %>% 
    unique() %>% 
    as.data.frame()
  colnames(fecal_16S_barcode_list) <- "barcode_list"
  
  write.table(fecal_16S_barcode_list, 
              file = "PRJEB22894_barcode_list.txt", 
              quote = F, row.names = F, col.names = F, sep = "\n")
  # 将sampleID和barcode分开
  fecal_16S_barcode_list_2 <- separate(fecal_16S_barcode_list, 
                                       barcode_list, 
                                       into = c("subject_id", "barcode"), 
                                       sep = "\\s") %>% 
    arrange(subject_id)
  # 查看不匹配的sample
  fecal_16S_barcode_list_2[!(fecal_16S_barcode_list_2$subject_id %in% 
                               sample_all$subject_id),]
}

# 整理样本信息mapping
if (T) {
  setwd("~/Wargo/EGAD00001003943/xmls/samples")

  filenames <- dir()
  
  #buccal_sample_name <- buccal@id %>% 
  #  as.character() %>% 
  #  str_match("(Wargo\\.[0-9]+)_.*") %>% .[,2] %>% unique() %>% 
  #  as.data.frame()
  
  #fecal_16S_sample_name <- fecal_16S@id %>% 
  #  as.character() %>% 
  #  str_match("(Wargo\\.[0-9]+)_.*") %>% .[,2] %>% unique() %>% 
  #  as.data.frame()
  #colnames(fecal_16S_sample_name) <- "subject_id"
  
  sample_all <- tibble()
  for (i in filenames) {
    meta <- xmlParse(file = i) %>% xmlRoot()
    meta_id <- xmlToList(meta[[1]][[4]]) %>% t() %>% as.tibble() 
    meta_data <- meta[[1]][[5]] %>% xmlToDataFrame() %>% t() %>% as.tibble()
    colnames(meta_data) <- meta_data[1,]
    meta_data <- meta_data[-1,]
    meta_extract <- bind_cols(meta_id, meta_data)
    sample_all <- bind_rows(sample_all, meta_extract)
  }
  colnames(sample_all)[1] <- "description"
  #sample_all <- sample_all %>% arrange(subject_id)
  #sample_all_fecal <- sample_all[sample_all$description == 
  #                                 "fecal sample 16S",]
  
  write.csv(sample_all, 
            file = "/home/yeguanhua/Wargo/Wargo_metadata.csv", 
            quote = F, row.names = F, na = "", append = F)
}
