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
Sampl <- read.csv("./Accession_list.csv")


###### Proc generation of the anatomies
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
}


# # write.csv(ROOT, "./GRANAR/anat.csv")
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
             kw = 0.00024,
             km = 3e-5,
             kAQP = 0.00043,
             kpl = 5.3e-12)

  wallthick(path = "MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml", 1.5)

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
fls <- list.files("./MECHA/Projects/GRANAR/out/M1v4/Root/res/")
fls <- fls[grepl(".txt", fls)]

K <- tibble(kr = NULL, kx = NULL, sampl_id = NULL, apo = NULL)
for (k in fls){
  M <- read_file(paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/res/",k))
  tmp_M <- strsplit(M, split="\n")[[1]]
  K_xyl_spec <- as.numeric(strsplit(tmp_M[15], " ")[[1]][5])
  kr_M <- as.numeric(strsplit(tmp_M[17], " ")[[1]][4])
  scenario <- round(parse_number(unlist(str_split(k,"_"))[3])/10)
  sampl_id <- parse_number(unlist(str_split(k,"_"))[4])
  K <- rbind(K, tibble(kr = kr_M, kx = K_xyl_spec, sampl_id = sampl_id, apo = scenario))
}

K = K %>%arrange(sampl_id, apo)

Sampl$kr1_new[Sampl$id %in% K$sampl_id] = K$kr[K$apo == 1]
Sampl$kr2_new[Sampl$id %in% K$sampl_id] = K$kr[K$apo == 2]
Sampl$kr3_new[Sampl$id %in% K$sampl_id] = K$kr[K$apo == 4]

Sampl$Kr1[Sampl$id %in% K$sampl_id] = K$kr[K$apo == 1]*2*pi*Sampl$radius[Sampl$id %in% K$sampl_id]
Sampl$Kr2[Sampl$id %in% K$sampl_id] = K$kr[K$apo == 2]*2*pi*Sampl$radius[Sampl$id %in% K$sampl_id]
Sampl$Kr3[Sampl$id %in% K$sampl_id] = K$kr[K$apo == 4]*2*pi*Sampl$radius[Sampl$id %in% K$sampl_id]

Acces = c("CONICO", "GORDO", "JALA", "NALTEL", "PALOME", "REVENT", "TABLON", "ZAPCHI")


Sampl%>%
  ggplot()+
  geom_point(aes(radius, kr1_new, colour = radius-r_stele), alpha = 0.5)+
  geom_point(aes(radius, kr1_new),shape = 1, colour = "red", size = 5,
             data = Sampl%>%filter(Accession %in% Acces))+
  geom_text(aes(x = radius+0.07*radius, y = kr1_new+kr1_new*0.04, label = Accession),
            angle = 45,
            check_overlap = TRUE,shape = 1, colour = "red", size = 5,
             data = Sampl%>%filter(Accession %in% Acces, !duplicated(Accession)))+
  geom_point(aes(radius, kr1_new),shape = 1, colour = "blue", size = 5,
             data = Sampl%>%filter(Accession %in% id_clust))+
  theme_dark()+
  labs(colour = "Cortex width [mm]")+
  ylab("Root Radial hydraulic conductivity [cm hPa-1 d-1] ")+
  xlab("Root Radius [mm]")+
  viridis::scale_colour_viridis()

Sampl%>%
  ggplot()+
  geom_point(aes(radius, Kr1, colour = radius-r_stele), alpha = 0.5)+
  geom_point(aes(radius, Kr1),shape = 1, colour = "red", size = 5,
             data = Sampl%>%filter(Accession %in% Acces))+
  geom_point(aes(radius, Kr1),shape = 1, colour = "blue", size = 5,
             data = Sampl%>%filter(Accession %in% id_clust))+
  theme_dark()+
  labs(colour = "Cortex width [mm]")+
  ylab("Root Radial hydraulic conductance [cm4 hPa-1 d-1] ")+
  xlab("Root Radius [mm]")+
  viridis::scale_colour_viridis()








