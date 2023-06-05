##########

# Heymans Adrien - june 2023
# Anatomical procedure - GRANAR
# hydraulic procedure - MECHA
# Script for the paper "RootDiversity" of McLaughlin et al. 2023

##########
print("loading dependencies")
library(tidyverse)
library(plyr)
library(deldir)
library(alphahull)
library(xml2)
library(sp)
library(viridis)
library(readxl)
`%!in%` <- compose(`!`, `%in%`)
# ----------------
print("loading GRANAR")
source("./GRANAR/R/granar.R")
source("./GRANAR/R/micro_hydro.R")
print("loading params")
params <- read_param_xml("./GRANAR/www/Zea_mays_CT.xml")
####
Sampl <- read.csv("./ConductivityData_AllFocalAccessions.csv")
# Unique set of hydraulic properties for the estimation of kr
Sampl$kw = 0.00024
Sampl$kAQP = 0.00043
Sampl$kpl = 5.3e-12
Sampl$thickness = 1.5
# uniform parameter for GRANAR
Sampl$radius = Sampl$root_radius
Sampl$RXA = pi*Sampl$radius^2
Sampl$r_stele = sqrt(Sampl$TSA/pi)
Sampl$TCA = Sampl$RXA-Sampl$TSA
Sampl$log_CW = log(Sampl$radius-Sampl$r_stele, exp(1))
Sampl$CF = exp(3.1091221+0.6718735*Sampl$log_CW)
Sampl$OneC = exp(Sampl$log_CW)/Sampl$CF
Sampl$nX = Sampl$NMV
Sampl$XMA_1 = Sampl$MVA/Sampl$NMV
Sampl$X_size= 2*sqrt(Sampl$XMA_1/pi)
Sampl$aerenchyma = Sampl$AA/Sampl$TCA


###### Proc generation of the anatomies
ROOT = NULL
for(i in 1:nrow(Sampl)){
  tmp_sampl <- Sampl[i,]  
  print(i)
  if(file.exists("./MECHA/cellsetdata/current_root.xml")){
    file.remove("./MECHA/cellsetdata/current_root.xml")
    file.remove("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml")
  }
  # Run GRANAR and change the parameter for the selected simulation
  sim <- run_granar(params, tmp_sampl)
  # Write the output
  file.copy("./MECHA/cellsetdata/current_root.xml", paste0("./MECHA/cellsetdata/root_",i,".xml"), overwrite = T)
  file.copy("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml",
            paste0( "./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_",i,".xml"), overwrite = T)
  tmp_root = sim$nodes
  ROOT = rbind(ROOT, tmp_root%>%mutate(Accession = tmp_sampl$Accession))
}


# write.csv(ROOT, "./GRANAR/anat.csv")
ROOT = read.csv("./GRANAR/anat.csv")

id_clust = unique(ROOT$Accession[grepl("cluster", ROOT$Accession)])
ROOT%>%
  mutate(Clust= ifelse(Accession %in% id_clust, T, F))%>%
  filter(Clust)%>%
  ggplot()+
  geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2))+
  theme_classic()+
  coord_fixed()+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  facet_wrap(~Accession)

# ggsave("~/GitHub/RootDiversity/img/anat_cluster.png", width=10, height=10)


# ### Proc: estimation of the radial hydraulic conductivities
fls <- list.files("./MECHA/cellsetdata/")
fls <- fls[grepl("root_", fls)]
for(j in fls){
  message("--------------")
  print(j)
  message("--------------")
  
  if(file.exists("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")){
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt")
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt")
    file.remove("./MECHA/cellsetdata/current_root.xml")
    file.remove("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml")
  }
  
  # Loading input files for the current estimation
  fc <- file.copy(paste0("./MECHA/cellsetdata/",j), "./MECHA/cellsetdata/current_root.xml", overwrite = T)
  if(fc == FALSE){next()}
  fc <- file.copy(paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_", parse_number(j), ".xml"),
                  paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml"), overwrite = T)
  if(fc == FALSE){next()}
  
  # MECHA input change
  id <- parse_number(j) 
  microhydro(path = "MECHA/Projects/GRANAR/in/Maize_hydraulics.xml",
             kw = Sampl$kw[parse_number(j)],
             km = Sampl$km[parse_number(j)], 
             kAQP = Sampl$kAQP[parse_number(j)],
             kpl = Sampl$kpl[parse_number(j)])
  wallthick(path = "MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml", Sampl$thickness[parse_number(j)])
  
  # Run MECHA - - - - - - - 
  system("python3 ./MECHA/MECHAv4_septa.py")
  message("python script has ended")
  
  # if works well, then:
  if(file.exists("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")){
    # Save output
    message ("success")
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_1,0_",id,".txt"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_2,1_",id,".txt"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_4,2_",id,".txt"), overwrite = T)
  }else{message ("fail and move to next simulation")}
  
}


# Read Mecha output
fls <- list.files("./MECHA/Projects/GRANAR/out/M1v4/Root/")
fls <- fls[grepl(".txt", fls)]

K <- tibble(kr = NULL, kx = NULL, sampl_id = NULL, apo = NULL)
for (k in fls){
  M <- read_file(paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/",k))
  tmp_M <- strsplit(M, split="\n")[[1]]
  K_xyl_spec <- as.numeric(strsplit(tmp_M[15], " ")[[1]][5])
  kr_M <- as.numeric(strsplit(tmp_M[17], " ")[[1]][4])
  scenario <- round(parse_number(unlist(str_split(k,"_"))[3])/10)
  sampl_id <- parse_number(unlist(str_split(k,"_"))[4])
  K <- rbind(K, tibble(kr = kr_M, kx = K_xyl_spec, sampl_id = sampl_id, apo = scenario))
}

K = K %>%arrange(sampl_id, apo)

Sampl$kr1_new = K$kr[K$apo == 1]
Sampl$kr2_new = K$kr[K$apo == 2]
Sampl$kr3_new = K$kr[K$apo == 4]

Sampl$Kr1 = K$kr[K$apo == 1]*2*pi*Sampl$root_radius
Sampl$Kr2 = K$kr[K$apo == 2]
Sampl$Kr3 = K$kr[K$apo == 4]

Sampl%>%
  ggplot()+
  geom_point(aes(kr_1_AA, kr1_new), colour = "blue")+ # Endodermal Casparian band
  geom_point(aes(kr_2, kr2_new), colour = "red")+ # Suberized endodermis
  geom_point(aes(kr_3, kr3_new), colour = "grey")+ # Suberized endodermis and Exodermal Casparian band
  geom_abline(slope = 1)+
  theme_classic()+
  ylab("kr Adrien")+
  xlab("kr Chloee")


  

  



