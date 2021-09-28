#Normalize metagenomic Hi-C data and detect spurious contacts using zero-inflated Negative Binominal regression frameworks
#Auther and maintainer: Yuxuan Du <yuxuandu@usc.edu>
#HiCzin R script depends on 'glmmTMB' package
library('glmmTMB')


HiCzin = function(contig_info_file , valid_contact_file , thres)
{
  sample_data = read.csv(valid_contact_file , header = F , sep = ',' )
  sample_data = as.data.frame(sample_data)
  colnames(sample_data) = c('index1' , 'index2' , 'contacts')
  contig_info = read.csv(contig_info_file , header = F , sep = ',' )
  contig_info = as.data.frame(contig_info)
  
  sample_data[ , 1] = sample_data[ , 1] + 1
  sample_data[ , 2] = sample_data[ , 2] + 1
  
  if(ncol(contig_info) == 4){
    colnames(contig_info) = c('contig_name' , 'site' , 'length' , 'coverage')
    contig_info[contig_info$site==0 , 2] = 1
    contig_info[contig_info$coverage==0 , 4] = min(contig_info[contig_info$coverage!=0 , 4])
    
    sample_len = rep(0 , nrow(sample_data))
    sample_site = rep(0 , nrow(sample_data))
    sample_cov = rep(0 , nrow(sample_data))
    
    for(i in 1:nrow(sample_data))
    {
      sample_site[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 2]) * 
                             as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 2]))
      
      sample_len[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 3]) * 
                            as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 3]))
      
      sample_cov[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 4]) * 
                            as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 4]))
    }
    
    sampleCon = as.numeric(sample_data[ , 3])
    
    mean_site = mean(sample_site)
    sd_site = sd(sample_site)
    mean_len = mean(sample_len)
    sd_len = sd(sample_len)
    mean_cov = mean(sample_cov)
    sd_cov = sd(sample_cov)
    
    sample_site = (sample_site-mean_site)/sd_site
    sample_len = (sample_len-mean_len)/sd_len
    sample_cov = (sample_cov-mean_cov)/sd_cov
    
    data_sample = cbind(sample_site , sample_len , sample_cov , sampleCon)
    data_sample = as.data.frame(data_sample)
    colnames(data_sample) = c('sample_site' , 'sample_len' , 'sample_cov' , 'sampleCon')
    
    tryCatch(
      {
        
        fit1 = glmmTMB(sampleCon~sample_site+sample_len+sample_cov, data = data_sample,
                       ziformula=~sample_site+sample_len+sample_cov , family=nbinom2)
        
      },
      error = function(e){
        message(e)
        message(paste("\nskip",  sep=" "))
      },
      warning = function(w){
        message(w)
        message(paste("\nskip",  sep=" "))
      }
    )
    coeff = as.numeric(fit1$fit$par)
    res_sample = sampleCon/exp(coeff[1] + coeff[2]*sample_site + coeff[3]*sample_len+ coeff[4]*sample_cov)
    index_nonzero = (res_sample > 0)
    res_sample_nonzero = res_sample[index_nonzero]
    perc = quantile(res_sample_nonzero , thres)
    result = c(coeff[1:4] , perc , mean_site , sd_site , mean_len , sd_len , mean_cov , sd_cov)
    return(result)
 }else{
    colnames(contig_info) = c('contig_name' , 'length' , 'coverage')
    
    sample_len = rep(0 , nrow(sample_data))
    sample_cov = rep(0 , nrow(sample_data))
    
    for(i in 1:nrow(sample_data))
    {
      
      sample_len[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 2]) * 
                            as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 2]))
      
      sample_cov[i] = log(as.numeric(contig_info[as.numeric(sample_data[i , 1]) , 3]) * 
                            as.numeric(contig_info[as.numeric(sample_data[i , 2]) , 3]))
    }
    
    sampleCon = as.numeric(sample_data[ , 3])
    
    mean_len = mean(sample_len)
    sd_len = sd(sample_len)
    mean_cov = mean(sample_cov)
    sd_cov = sd(sample_cov)
    
    sample_len = (sample_len-mean_len)/sd_len
    sample_cov = (sample_cov-mean_cov)/sd_cov
    
    data_sample = cbind(sample_len , sample_cov , sampleCon)
    data_sample = as.data.frame(data_sample)
    colnames(data_sample) = c('sample_len' , 'sample_cov' , 'sampleCon')
    
    tryCatch(
      {
        
        fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                       ziformula=~sample_len+sample_cov , family=nbinom2)
        
      },
      error = function(e){
        message(e)
        message(paste("\nskip",  sep=" "))
      },
      warning = function(w){
        message(w)
        message(paste("\nskip",  sep=" "))
      }
    )
    
    coeff = as.numeric(fit1$fit$par)
    res_sample = sampleCon/exp(coeff[1] + coeff[2]*sample_len + coeff[3]*sample_cov)
    index_nonzero = (res_sample > 0)
    res_sample_nonzero = res_sample[index_nonzero]
    perc = quantile(res_sample_nonzero , thres)
    result = c(coeff[1:3] , perc , mean_len , sd_len , mean_cov , sd_cov)
    return(result)
  }
}













