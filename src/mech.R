


### loading dependencies
library(tidyverse)
library(plyr)
library(deldir)
library(alphahull)
library(xml2)
library(sp)
library(viridis)
library(readxl)
`%!in%` <- compose(`!`, `%in%`)

#### loading GRANAR 
source("./GRANAR/R/granar.R")
source("./GRANAR/R/micro_hydro.R")
### loading mock parameter file
params <- read_param_xml("./GRANAR/www/Zea_mays_CT.xml")
#### loading input from anatomical trait dataset
Raw_data <- readxl::read_excel("./www/Root no CR6_update_2.xlsx")

Raw_data[Raw_data$Treatment=='Drought',2]<-'sheltered'
Raw_data[Raw_data$Treatment=='Well watered',2]<-'non-sheltered'

Sampl = Raw_data%>%
  mutate(RXA = `Root area`,
         TSA = `Stele area`,
         TCA = `Cortex area`,
         AA = `Aerenchyma area`,
         aerenchyma= `Aerenchyma percent`/100,
         radius = sqrt(RXA/pi),
         r_stele = sqrt(TSA/pi),
         MXA = `Total Metaxylem area`,
         nX = `Metaxylem number`,
         X_size = 2*sqrt((MXA/nX)/pi),
         CF = `Cortical file number`,
         OneC = (radius-r_stele)/(CF+2),
         OC = 2*sqrt(`Cortical cell size`/pi),
         ratio = (2+0.07456*r_stele*1000)/nX,
         nPX = nX*ratio,
         PXA_1 = 1000^2*(sqrt(radius/35)/10)^2,
         k_protxyl_s = PXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
         kx_unM = k_protxyl_s*nPX*200/1E4, # kx when only the proto xylem have their cell wall lignified 
         LMXA = MXA,
         LMXA_1 = LMXA*1000^2/nX,
         k_Mxyl_s = LMXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
         kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM)

Sampl%>%
  ggplot()+geom_histogram(aes(RXA-TCA), fill = 'blue', alpha = 0.5, bins = 50)+geom_histogram(aes(TSA), fill = 'red', alpha = 0.5, bins = 50)+
  theme_classic()+xlab('TCA [mm2]')

av_dat = Sampl%>%
  dplyr::group_by(Genotype, Treatment, Roottype)%>%
  dplyr::summarise(m_CF = mean(CF, na.rm = T),
                   m_OC = mean(OC, na.rm = T),
                   m_aerenchyma = mean(aerenchyma, na.rm = T), .groups = "drop")

Sampl = left_join(Sampl, av_dat , by = c("Genotype", "Treatment", "Roottype"))

Sampl = Sampl%>%
  mutate(CF = ifelse((is.na(CF) | CF <= 4), m_CF,CF),
         OneC = (radius-r_stele)/CF,
         OC = ifelse((is.na(OC) | OC > 0.2), m_OC,OC),
         aerenchyma = ifelse(is.na(aerenchyma), m_aerenchyma,aerenchyma))%>%
  arrange(radius)

Sampl$id = 1:nrow(Sampl)


# ### Proc: estimation of the radial hydraulic conductivities
fls <- list.files("./MECHA/cellsetdata/")
fls <- fls[grepl("root_", fls)]

fl_done = list.files("./MECHA/Projects/GRANAR/out/M1v4/Root/")
fl_done = fl_done[grepl("Macro_prop_1,0_", fl_done)]
ids = parse_number(unlist(str_split(fl_done,"1,0")))
ids <- ids[!is.na(ids)]
fl_done = paste0("root_", ids,".xml")
fls = fls[fls %!in% fl_done]


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