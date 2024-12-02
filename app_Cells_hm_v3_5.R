#### app_Cells_hm_v3.R : added annotation selection, change cluster order etc 
## 3.3 : corrected gene log heatmap to normalize first by cpm the gene expression by metacell
## 3.4 adds PlasmaCells, and support for tglkmeans clusering
## 3.5 adds support for personalized barplots
# Version of the script : 30 Nov 2024

# THOMAS : 
googlepath="G:/Mon Drive"

# JEROME :
googlepath="~/Google Drive"

#MONIKA:
googlepath="G:/.shortcut-targets-by-id/14-qt66Lhuv7qgN9Rf0BiBSPPp0_bfPCG"

setwd(paste0(googlepath,"/UG_metacells/McExplorer_cell_lineages"))

##################################
install.packages("BiocManager")
install.packages("devtools")
install.packages(c("shiny","ggplot2","magick","RColorBrewer","dplyr","seriation"))
install.packages(c("raster","dendsort","matrixStats","datasets","circlize","heatmaply","plotly","scales","Matrix")             
install.packages("ComplexHeatmap")
install.packages("shinyWidgets")
install.packages(c("patchwork","tidyverse"))
##################################



library(shiny)
library(ggplot2)
library(magick)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(raster)
library(seriation)
library(dendsort)
library(matrixStats)
library(datasets)
library(circlize)
library(heatmaply)
library(plotly)
library(scales)
#library(pryr)  # Nico l'a pas, a quoi il sert ?
library(Matrix)
#library(Rfast)
library(shinyWidgets)
library(patchwork)
library(tidyverse)

# Load the data depending of the cell type to explore:
################################################################################
palette_en_cours<-rep("black",50)

############
# FOR Macs:
object="MacsfromMNP"
load("genelist_16nov22_Macs.rd") # genelist 150 Modules found on Enriched Macs
load(paste0(googlepath,"/UG_metacells/McExplorer_cell_lineages/MacsTot_z_mymod_150_metadata_goodlog_v3.rd"))
if(file.exists(paste0(object,"_select2_dgCmc_counts.rd"))){
  load(paste0(object,"_select2_dgCmc_counts.rd")) # dgcmtxcounts (has the metacells/gene counts)
}else{ load(paste0(object,"_dgCmc_counts.rd")) }# dgcmtxcounts (has the metacells/gene counts)
# Pour avoir l'ordre des modules comme dans les Enrichis :
load(paste0(googlepath,"/UG_metacells/Figures/ht_Enrich_Macs.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Total_Macs.rd"))


###########
# FOR DCs:
object="DCfromMNP"
load("genelist_16nov22_DC.rd")
load(paste0(googlepath,"/UG_metacells/McExplorer_cell_lineages/DCTot_z_mymod_150_metadata_goodlog_v3.rd"))
if(file.exists(paste0(object,"_select2_dgCmc_counts.rd"))){
  load(paste0(object,"_select2_dgCmc_counts.rd")) # dgcmtxcounts (has the metacells/gene counts)
} else{ load(paste0(object,"_dgCmc_counts.rd")) }# dgcmtxcounts (has the metacells/gene counts)
# Pour avoir l'ordre des modules comme dans les Enrichis :
load(paste0(googlepath,"/UG_metacells/Figures/ht_Enrich_DC.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Total_DC.rd"))


###########
# For Fibro:
object="Totdec_StrGli_Fibro"
load(paste0("SubsetsAppObj/genelist_150m_",object,".rd"))
load(paste0("SubsetsAppObj/",object,"_z_150m_metadata.rd"))
load(paste0("SubsetsAppObj/",object,"_dgCmc_counts.rd"))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Total_Fibro.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))

###########
# For B cells:
object="Totdec_B"
load(paste0("SubsetsAppObj/genelist_150m_",object,".rd"))
load(paste0("SubsetsAppObj/",object,"_z_150m_metadata.rd"))
load(paste0("SubsetsAppObj/",object,"_dgCmc_counts.rd"))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Totdec_B.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))

###########
# For Tcells:
object="Totdec_Tcellsmini"
load(paste0("SubsetsAppObj/genelist_150m_",object,".rd"))
load(paste0("SubsetsAppObj/",object,"_z_150m_metadata.rd"))
load(paste0("SubsetsAppObj/",object,"_dgCmc_counts.rd"))
#load(paste0(googlepath,"/UG_metacells/Figures/ht_",object,".rd"))
# Corrected HM
load("SubsetsAppObj/babo_zmodMM_Totdec_Tcellsmini.rd") # Contient les modules corrig?s
# babo_zmodbis contains corrected modules in MyModXsca 
load(paste0(googlepath,"/UG_metacells/Figures/ht_Totdec_Tcellsmini_corrected.rd"))
#palette_en_cours<-c('T_1 ILC3-like'='#874037','T_2 Innate-like'='#A3672C','T_3 Innate-like CD8'='#C3A34B','T_4 CD8 cytotox'='#D6D893','T_5 CD8 Trm'='#4F88B9','T_6 CD4 Trm'='#74BBCD','T_7 Trm conta macs?'='#a4dceb',
#                   'T_8 naive/TCM'='#f77774','T_9 Treg'='#2ca02c','T_10 TFH_CXCL13'='#d29fe3','T_11 Prolif T'='#e151e3','T_12 B conta'='#BCBCBB' )

###########
# For Plasma cells:
object="Totdec_PlasmaCells"
load(paste0("SubsetsAppObj/genelist_150m_",object,".rd"))
load(paste0("SubsetsAppObj/",object,"_z_150m_metadata.rd"))
load(paste0("SubsetsAppObj/",object,"_dgCmc_counts.rd"))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Totdec_PlasmaCells.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))

###########
# For EndoC (with Lymphatics):
object="Totdec_StrGli_EndoC"
load(paste0("SubsetsAppObj/genelist_150m_",object,".rd"))
load(paste0("SubsetsAppObj/",object,"_z_150m_metadata.rd"))
load(paste0("SubsetsAppObj/",object,"_dgCmc_counts.rd"))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Totdec_StrGli_EndoC.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))

###########
# For Glial Pericytes SmMuscle (with Lymphatics):
object="Totdec_StrGli_GliPerSM"
load(paste0("SubsetsAppObj/genelist_150m_",object,".rd"))
load(paste0("SubsetsAppObj/",object,"_z_150m_metadata.rd"))
load(paste0("SubsetsAppObj/",object,"_dgCmc_counts.rd"))
load(paste0(googlepath,"/UG_metacells/Figures/ht_Totdec_StrGli_GliPerSM.rd"))
row_mod_ord<- as.vector(unlist(row_order(ht)))


#---------------------------------------#
rm(km_k40)
if (file.exists(paste0("./SubsetsAppObj/km_k40_",object,".rd") )) {
  load(paste0("./SubsetsAppObj/km_k40_",object,".rd"))
} else {
  print("reordering file does not exist")
}
#---------------------------------------#



colnames(genelist)<-c("Mods","Gens")
hm_zmod<-as.matrix(metadf_z[,c(154:303)])
df_zmod2<-as.data.frame(metadf_z[,c(2:3,154:303)])
colnames(df_zmod2)<-sub("([a-zA-Z_]+)(\\d+)([a-zA-Z]+)", "Module \\2", colnames(df_zmod2))  
#colnames(hm_zmod)<-gsub("scaled","",colnames(hm_zmod))

if(object=="Totdec_Tcellsmini"){
  # For T cells, using corrected modules as hm_zmod:
  hm_zmod<-as.matrix(babo_zmodbis[grep("^MyMod[1-9]\\d*sca$", colnames(babo_zmodbis), value = TRUE)] )
  colnames(hm_zmod)<-sub("([a-zA-Z]+)(\\d+)([a-zA-Z]+)", "My_Mod_\\2\\3led", colnames(hm_zmod))
  df_zmod2<-as.data.frame(babo_zmodbis[,c(2:3,304:306)])
  df_zmod2<-cbind(df_zmod2,hm_zmod)
  colnames(df_zmod2)[6:155]<-sub("([a-zA-Z_]+)(\\d+)([a-zA-Z]+)", "Module \\2", colnames(df_zmod2)[6:155])  #Oblig? de mettre chiffre pour pas qu'il renom les V1pct etc
  
}

preloco_nc<-as.matrix(log(1+dgcmtxcounts))# DU coup les matrices de genes sont pas corrig?es, on peut pas
allsum<-Matrix::colSums(preloco_nc) 
row_mod_ord<- as.vector(unlist(row_order(ht)))
data.cpm<-1e+06*sweep(dgcmtxcounts,2,colSums(dgcmtxcounts),`/`)  # if you don't want another package
# l.data.cpm<-log(1+data.cpm)
# library(edgeR)
# eldata.cpm<-cpm(dgcmtxcounts, normalized.lib.sizes=TRUE, log=T,prior.count = 1)
# 
# preloco<-as.matrix(l.data.cpm)
# preloco<-as.matrix(eldata.cpm)
preloco<-as.matrix(log2(1+(data.cpm/10)))  # Divide by 10 the CPM for more readability

# as.matrix(eldata.cpm)[45:49,1:12]-as.matrix(eldata.cpm)[5,1] #is almost the same as:
# log2(1+(data.cpm/10))[45:49,1:12]
# 
# mtX<-as.matrix(eldata.cpm)[5:9,1:6]-as.matrix(eldata.cpm)[5,1]
# mtY<-as.matrix(log(1+dgcmtxcounts))[5:9,1:6]
# sweep(mtX,1,mtY,FUN="/")


###----------------------------------------------------------### 
#### Function to transform list of genes to vector readable #### 
#### by the program, with exclusion of unknown genes        ####
###----------------------------------------------------------### 
gene_input2vect<-function(genelist,dgcmtx){
  gnshow<-unlist(strsplit(genelist,split=","))
  gnshow<-gsub(" ","",gnshow)
  notfound<-setdiff(gnshow,rownames(dgcmtx))
  print( paste0("[",paste0(notfound,collapse =","),"] not found") )
  return(setdiff(gnshow,notfound) )
}
####################################################

#ifexist(){genelist=gnlscenics}

findgene <- function(x,y){
  for(i in 1:length(x)){
    if(y%in%unlist(x[i])){return(i)}
  }
}
# List for the MAcrophage figure: gnlist<-"CXCL8,SERPINB2,MFSD2A,F3,IL1A,TNF,IL7R,IL3RA,IDO1,INHBA,CD274,P2RX7,SLAMF1,CRIM1,SLC39A8,SLC7A11,ADAM19,SPP1,DNAJB4,HSPA6,FCGR2B,CA2,SCD,BNIP3,RRAD,SLC2A1,IFI44L,CXCL10,CXCL9,IFIT2,OAS2,CXCL11,IL23A,CD40,ICAM4,SCN1B,SLAMF7,CCL20,CXCL3,CXCL2,S100A8,S100A9,VCAN,FCN1,RIPOR2,S100A12,APOBEC3A,SELL,FCAR,ALOX5AP,RETN,MCEMP1,OSM,ATP5PB,MRPS21,RNF7,NDUFB5,ATP5ME,C1QA,VSIG4,DNASE1L3,DAB2,PDK4,SLCO2B1,CMKLR1,CD209,FUCA1,MERTK,SELENOP,FOLR2,IGF1,AXL,LYVE1"
gnlist<-"CD7,CD2,CD3D,STMN1,PCNA,TNFRSF4,TNFRSF18,CTLA4,IL2RA,FOXP3,CCR7,SELL,LEF1,CST7,GZMA,GZMK,CD8A,CD8B,CD69,ID2,ANKRD28,GPR171,TNFAIP3,ANXA1,IL7R,TRDC,TYROBP,FCER1G,PTGDR,KRT81,KRT86,KIT,PCDH9,ALDOC,LINC00299,CMC1, XCL2, IL2RB, KLRD1,KLRF1, CLIC3, MCTP2,MZB1,TNFRSF17,SEC11C,DERL3,XBP1,IGHA2,IGHM,IGHG1,IGHG3,BANK1, VPREB3, CD24, ARHGAP24, FCRLA, RALGPS2, TNFRSF13C, SPIB,MS4A1, CD79B, CD19,IGHD, FCER2,TCL1A, CD72,CD27,CLECL1,TNFRSF13B,FCGR3A,CD14,CD68, SLC40A1, STAB1, CSF1R, MS4A4A, SLCO2B1, MAFB,MS4A7,C1QA,C1QB,C1QC,IL1RN,S100A8,S100A9,TREM1,KYNU,C1orf54,CLEC9A,CADM1,CPNE3,XCR1,IDO1,CLEC10A,FCER1A,CD1C,CD1D,GPR157,LAMP3,DAPP1,FSCN1,CCL19,CCL22,EBI3,GZMB,TCF4,LILRA4,IRF7,CLEC4C,IL1RL1,CPA3,MS4A2,TPSB2,TPSAB1,ADCYAP1,GATA2,HPGDS,HPGD,VWF, RAMP3, NPDC1, JAM2, PLVAP, NOTCH4, HSPG2, ESAM, CYYR1, CD93, ICAM2, S1PR1, RAMP2,CD36, CA4, TMEM88, FLT1,ACKR1, SELP,DUSP23,PROX1,MMRN1,CCL21,IGF2,LYVE1,RGS5,NDUFA4L2,C11orf96,ACTG2,MYH11,CXCL14,ADH1B,CTSK,MMP2,LUM,PTGDS,TCF21,ADAMDEC1,CCL13,CCL8,ADAM28,HAPLN1, ABCA8,CFD,THY1,CHI3L1,BGN,PDPN,TNFRSF12A,MPZ,NRXN1,SCN7A,SEMA3B,SOX2,MYOT,GFRA3,TUBB2B,GPM6B,PLP1,XKR4,ALDH1A1,CRYAB,NTM,ANK3,LGI4,S100B"
ann_df<-data.frame(metacell=rownames(metadf_z),annotation="undefined")
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  navbarPage("MCexplorer_HM v3.5 IBD",
             tabPanel("z-scored HM",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Heatmap of the Z-scaled module, displayed as mean in metacells"),
                                 verbatimTextOutput('prefoo')
                               )
                        ),
                        mainPanel(
                          plotOutput("hmzall",brush = "ht_brush", click = "ht_click", height = 500, width = 900),
                          tags$hr()
                        )
                      )
                      
             ),
             tabPanel("z-scored Heatmaply",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Heatmap of the Z-scaled modules enrichment by metacells"),
                                 textAreaInput("clsordhm", h4("reorder clusters (n?)"),height ='200%',width='100%',
                                               value = paste0(c(1:40),collapse=",")),
                                 selectInput("var", label = "Choose a module to display", choices=c(1:150),selected=1),
                                 textInput("geneask", h4("Looking for a gene ?"),placeholder = "no gene selected"),
                                 verbatimTextOutput("vgen"),
                                 # sliderInput(inputId = "HMwidz2", label = "Width", value=1000, min=0, max=1600,step=100),
                                 # sliderInput(inputId = "HMheign2", label = "Height", value=800, min=0, max=1600,step=100),
                                 sliderInput(inputId = "hmgnprop", label = "Gene heatmap size", value=0.4, min=0, max=1,step=0.05),
                                 sliderInput(inputId = "hmgnlbsz", label = "HM y-labels size", value=0.75, min=0, max=2,step=0.05),
                                 textAreaInput("geneshow", h4("Genes to plot"),height ='200%',width='100%',
                                               value = gnlist),
                                 materialSwitch("interact_on", value=F,label = HTML(paste0(
                                   "Activate interaction\n <span style='color: red;font-weight: bold;'> /!\\ It will be slower</span> "
                                   )), status = "danger"),
                                 materialSwitch("reorder_on", value=F,label = HTML(paste0(
                                   "Reorder from tglk file"
                                 )), status = "danger"),
                                 verbatimTextOutput("reor1"),
                                 # radioButtons("radio_inter", label = HTML(paste0(
                                 #      "Activate interaction\n <span style='color: red;font-weight: bold;'>/!\\ It will be slower</span> ", list("Nope", "yes"), 
                                 #              selected = "Nope"),
                                 actionButton("click", "Run Heatmap"),
                                 tags$hr(),
                                 verbatimTextOutput("notf_gen"),
                                 textInput("clusterschi", h4("Chi? test"),value = ""),
                                 actionButton("click_run_chi", "Run chi2 test"),
                                 verbatimTextOutput("chirep"),
                               )
                        ),
                        mainPanel(
                          textOutput("text"),
                          verbatimTextOutput("verb"),
                          tags$hr(),
                          #plotlyOutput("hmzall2", height = 1500, width = 1000),
                          conditionalPanel(
                            condition = "input.interact_on == false",
                            plotOutput("hmzall2_gg", height = 1500, width = 1000)
                          ),
                          conditionalPanel(
                            condition = "input.interact_on == true",
                            plotlyOutput("hmzall2_inter", height = 1500, width = 1000)
                          ),
                          #tags$hr(),
                          ##    verbatimTextOutput("test"),
                          #plotlyOutput("hmgene2", height = 400, width = 900),
                          textInput("topmod", h4("Top modules (>0.8) for the cluster:"),value = ""),
                          verbatimTextOutput("topmodout"),
                        )
                      )   
             ),
             tabPanel("Modscores by clusters",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Heatmap of the Z-scaled modules enrichment by metacells"),
                                 selectInput("var2", label = "Choose a module to display", choices=c(1:150),selected=1),
                               )
                        ),
                        mainPanel(
                          textOutput("text2"),
                          verbatimTextOutput("verb2"),
                          plotOutput("mod_clust_bar", height = 500, width = 900),
                        )
                      )
             ),
             tabPanel("Mc annotation",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Annotation"),
                      fluidRow(column(4,
                                      wellPanel(
                                        h3("Cluster selection"),
                                        textAreaInput("ann_clu",label= NULL,height ='200%',width='100%',placeholder="Cluster numbers, comma-separated")
                                      )),
                               column(4,
                                      wellPanel(
                                        h3("Add metacells"),
                                        fluidRow(
                                          column(2, HTML('<h4><b>Mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_add_mc1",label= NULL,placeholder="XXX"),offset = 0, style='padding:0px;'),
                                          column(3, HTML('<h4><b>to mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_add_mc2",label= NULL,placeholder="XXX"), offset = 0,style='padding:0px;')
                                        )
                                      )),
                               column(4,
                                      wellPanel(
                                        h3("Remove metacells"),
                                        fluidRow(
                                          column(2, HTML('<h4><b>Mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_rem_mc1",label= NULL,placeholder="XXX"),offset = 0, style='padding:0px;'),
                                          column(3, HTML('<h4><b>to mc</b></h4>'),offset = 0, style='padding-right:0px;'),
                                          column(3, textInput(inputId = "ann_rem_mc2",label= NULL,placeholder="XXX"), offset = 0,style='padding:0px;')
                                        )
                                      ))
                      ),
                      #br(),
                      #actionButton("click", "get metacell list"),
                      tags$hr(),
                      fluidRow(
                        column(8,
                               wellPanel(
                                 h3("Selected metacells to annotate"),
                                 verbatimTextOutput("out_mcs"),
                               )),
                        column(4,
                               wellPanel(
                                 h4("Name of the group"),
                                 textInput(inputId = "grp_name",label= NULL,value=""),
                                 actionButton("btn_save_name", "save name"),
                               ))
                      ),
                      fluidRow(
                        column(8,
                               wellPanel(
                                 h3("Table of proportions"),
                                 tableOutput("out_prop_cellannot")
                               )),
                        column(4,
                               wellPanel(
                                 textInput("file_name", "File Name (without extension):"),
                                 downloadButton("save_csv", "Save CSV")
                               ))),
                      fileInput("annot_csv", "Or choose an annotation file:"),
                      actionButton("btn_load_name", "load name")
             ),
             tabPanel("Modscores by group",
                      tags$head(tags$style(HTML('* {font-family: "Roboto"};'))),
                      titlePanel("Module selection"),
                      fluidRow(
                        column(3, 
                               wellPanel(
                                 h4("Boxplot of the Z-scaled modules enrichment or gene expression by cell category"),
                                 # selectInput("var3", label = "Choose a module to display", choices=c(1:150),selected=1),
                                 selectInput("name_fromlist", "Select data to plot \n(personal gene signatures are un-corrected for Tcells)", choices = c("Gene (log CPM) or Gene signature",colnames(df_zmod2))),
                                 uiOutput("reactive_txt"),
                                 textInput("palette_1", "personnalized palette",value = "" ),
                                 textInput("bx_order", "X axis order", value=""),
                                 switchInput("toggle_legend", "Show Legend", value = TRUE),
                                 sliderInput(inputId = "bx_hei", label = "Height", value=600, min=0, max=1200,step=50),
                                 sliderInput(inputId = "bx_wid", label = "Width", value=550, min=0, max=1200,step=50),
                                 sliderInput(inputId = "bx_siz", label = "Size", value=1.3, min=0, max=3,step=0.1)
                                 # switchInput("toggle_code", "Show Code", value = TRUE)
                                 )
                               
                        ),
                        mainPanel(
                          textOutput("text3"),
                          verbatimTextOutput("verb3"),
                          plotOutput("mod_group_bar", height = 500, width = 900),
                        )
                      )
             )
             
  ))

hm_env = new.env()

server <- function(input, output, session) {
  ## Search and provide genes from modules:
  output$text <- renderText({paste0(genelist$Mod[as.numeric(input$var)]," genes: ")})
  output$verb <- renderText({genelist$Gens[as.numeric(input$var)] })
  output$text2 <- renderText({paste0(genelist$Mod[as.numeric(input$var2)]," genes: ")})
  output$verb2 <- renderText({genelist$Gens[as.numeric(input$var2)] })
  output$text3 <- renderText({paste0(genelist$Mod[as.numeric(sub("^Module ","",input$name_fromlist))]," genes: ")})
  output$verb3 <- renderText({genelist$Gens[as.numeric(sub("^Module ","",input$name_fromlist))] })
  output$vgen <- renderText({
    if(input$geneask%in%unlist(strsplit(genelist$Gens,split=","))){
      print(paste0("Module ",findgene(strsplit(genelist$Gens,split=","),input$geneask)))
    }else{print("Not found")}
  })
  bxhei<-function(){return(input$bx_hei)}
  bxwid<-function(){return(input$bx_wid)}
  bxsiz<-function(){return(input$bx_siz)}
  
  output$hmzall <- renderPlot({
    print("lezgo")
    ##Compute the heatmap figure (Macs/DC have different code as we already have a fixed order to match our previous analyses)
    set.seed(42)
    validate(need(exists("ht")==T,"HT OBJECT NOT FOUND / Please precompute heatmap \n\n(it allows a common ground of matrix order for all configs \nas setting seed doesn't always work)"))
    hm_env$ht = ht 
    ComplexHeatmap::draw(hm_env$ht)  # semblerait plus long de l'afficher en appelant draw ??? mais que hors shiny, dedans c'est pareil ?
    #print(hm_env$ht)
    print("draw")
    
    
    metacells_kclust<- column_order(hm_env$ht) 
    hm_env$metacells_kclust<-metacells_kclust
    
    for (i in 1:length(metacells_kclust)){
      if (i == 1) {
        clu <- t(t(rownames(hm_zmod)[metacells_kclust[[i]]]))
        out <- cbind(clu, paste("cluster", i, sep=""))
        colnames(out) <- c("Metacell", "Cluster")
      } else {
        clu <- t(t(rownames(hm_zmod)[metacells_kclust[[i]]]))
        clu <- cbind(clu, paste("cluster", i, sep=""))
        out <- rbind(out, clu)
      }
    }
    rm(clu)
    out<-as.data.frame(out)
    hm_env$out<-out
    
    # Automatic Table to display on the shiny app (Later do a tool to select a k cluster (can interactive give which cluster?) and it shows the list the genes in Module positives) :
    mod_clu<-c()
    mcs_clu<-list()
    nmod_clu<-list()
    for(i in 1:40){
      clust=paste0("cluster",i)
      mcs_clu[[i]]<-out$Metacell[which(out$Cluster%in%clust)]
      mod_clu[[i]]<-if(length(mcs_clu[[i]])==1){colnames(hm_zmod)[which(hm_zmod[mcs_clu[[i]],]>0.8)]
      }else{colnames(hm_zmod)[which(colMeans2(as.matrix(hm_zmod[mcs_clu[[i]],]))>0.8)]}
      nmod_clu[[i]]<-length(mod_clu[[i]])
    }
    hm_env$mcs_clu<-mcs_clu
    hm_env$mod_clu<-mod_clu
    })


  output$notf_gen = renderText({
    gnshow<-unlist(strsplit(input$geneshow,split=","))
    gnshow<-gsub(" ","",gnshow)
    notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
    print( paste0("[",paste0(notfound,collapse =","),"] not found") )
    ## ## gc()
  })
  
  reor_switch_state <- reactive({
    input$reorder_on
  })
  int_switch_state <- reactive({
    input$interact_on
  })
  
  output$hmzall2_gg <- renderPlot({
    
    input$click
    req(input$click) #to prevent print at first lauch
    isolate({
      print("plotting ggplot HM")
      clorhm<-as.numeric(unlist(strsplit(input$clsordhm,split=",")))
      htmxord<-t(hm_zmod)[ as.vector(unlist(row_order(hm_env$ht))) , as.vector(unlist(column_order(hm_env$ht)[clorhm])) ] 
      #htmxord<-t(hm_zmod)[ as.vector(unlist(row_order(hm_env$ht))) , as.vector(unlist(column_order(hm_env$ht))) ] 
      gnshow<-unlist(strsplit(input$geneshow,split=","))                                                                          
      gnshow<-gsub(" ","",gnshow)                                                                                                 
      notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
      loco<-preloco[setdiff(gnshow,notfound),colnames(htmxord)]
      #If no interaction selected, print the ggplots, else, print the heatmaply
      
      switch_reor <- reor_switch_state()
      switch_int <- int_switch_state()
      
     #v3.3#  if (input$interact_on==FALSE) {
      if (switch_int==FALSE) {
        hm_env$hmgg1<-htmxord %>% 
          as.data.frame() %>%
          rownames_to_column("modules") %>%
          pivot_longer(-c(modules), names_to = "metacells", values_to = "counts")
        # format(object.size(htmxord),units="auto")  2.1Mb    format(object.size(hmgg1),units="auto")  6.1Mb 
        # hmgg1 est 3 fois plus gros que htmxord
        
        if (switch_reor==TRUE) {
          lim_ord=km_k40$cluster$id[order(km_k40$cluster$clust)]}else{
            lim_ord=colnames(htmxord)
          }
        
        hm_env$ggpt1<-ggplot(hm_env$hmgg1,aes(x=metacells, y=modules, fill=counts)) + 
          geom_raster() + scale_x_discrete(limits=lim_ord)+ 
          scale_y_discrete(limits=rev(rownames(htmxord)))+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y = element_text(size = rel(input$hmgnlbsz) ) )+
          scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish)
        
        hm_env$hmgg2<-loco%>% 
          as.data.frame() %>%
          rownames_to_column("genes") %>%
          pivot_longer(-c(genes), names_to = "metacells", values_to = "counts")
        # format(object.size(loco),units="auto")  2.6Mb    format(object.size(hmgg2)),units="auto")  7.4Mb 
        # hmgg2 est 3 fois plus gros que loco
        # if(length(setdiff(gnshow,notfound))<10){fntsiz<-1
        # }else if(length(setdiff(gnshow,notfound))<30){fntsiz<-0.8
        # }else{fntsiz<-0.6}
        hm_env$ggpt2<-ggplot(hm_env$hmgg2,aes(x=metacells, y=genes, fill=counts)) + 
          geom_raster() + scale_x_discrete(limits=lim_ord)+ 
          scale_y_discrete(limits=rev(rownames(loco)))+
          scale_fill_gradientn(colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis"))+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y = element_text(size = rel(input$hmgnlbsz)) )
       print(hm_env$ggpt1/hm_env$ggpt2+plot_layout(heights = c((1-input$hmgnprop),input$hmgnprop)) )
      }
        
    })
    
  })
      
  output$hmzall2_inter  = renderPlotly({
    
    input$click
    req(input$click) #to prevent print at first lauch
    isolate({
      
      print("plotting heatmapply HM")
      clorhm<-as.numeric(unlist(strsplit(input$clsordhm,split=",")))
      htmxord<-t(hm_zmod)[ as.vector(unlist(row_order(hm_env$ht))) , as.vector(unlist(column_order(hm_env$ht)[clorhm])) ] 
      #htmxord<-t(hm_zmod)[ as.vector(unlist(row_order(hm_env$ht))) , as.vector(unlist(column_order(hm_env$ht))) ] 
      gnshow<-unlist(strsplit(input$geneshow,split=","))                                                                          
      gnshow<-gsub(" ","",gnshow)                                                                                                 
      notfound<-setdiff(gnshow,rownames(dgcmtxcounts))
      loco<-preloco[setdiff(gnshow,notfound),colnames(htmxord)]
      
     ##v3.3#if (input$interact_on==TRUE) {
      switch_reor <- reor_switch_state()
      switch_int <- int_switch_state()
      
      if (switch_int==TRUE) {
        if (switch_reor==TRUE) {
          lim_ord=km_k40$cluster$id[order(km_k40$cluster$clust)]}else{
            lim_ord=colnames(htmxord)
          }
        
      ## Prepare the same heatmap as in 1st window but not cut in clusters, and interactive (plotly)
      hm_env$htpt1<-heatmaply(htmxord[,lim_ord], dendrogram = "none",showticklabels = c(F,T),limits=range(-2,2), fontsize_row = rel(input$hmgnlbsz),colors = "RdBu",Rowv=F,Colv=F,legendgroup="1st",showlegend = T,coloraxis = 'coloraxis',
                              scale_fill_gradient_fun =scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits=c(-2, 2), oob=squish) )


      ##  Second Heatmap displaying the selected genes (below the module heatmap)

      ###loco<-as.matrix(log(1+dgcmtxcounts[setdiff(gnshow,notfound),colnames(htmxord)]))
      # old code to have a line showing the PBMC vs Tissue distribution, might be usefull later to show other infos ?
      # pt1.5<-heatmaply(t(mtxtot), dendrogram = "none",limits=range(0,1),colors = "RdBu",Rowv=F,Colv=F,legendgroup="1st",showlegend = T,coloraxis = 'coloraxis',
      #                  scale_fill_gradient_fun =scale_fill_gradient2(low = "blue",mid="white", high = "red", limits=c(0, 1),midpoint = 0.5) )
      
      
      # if(length(setdiff(gnshow,notfound))<10){fntsiz<-8
      # }else if(length(setdiff(gnshow,notfound))<30){fntsiz<-7
      # }else{fntsiz<-6}
      hm_env$htpt2<-heatmaply(loco[,lim_ord], fontsize_row = rel(input$hmgnlbsz),fontsize_col = 7,showticklabels = c(F,T),dendrogram = "none",Rowv=F,Colv=F,legendgroup="2nd",coloraxis = 'coloraxis2',scale_fill_gradient_fun =scale_fill_gradientn(
        colors = viridis(n = 256, alpha = 1, begin = 0,end = 1, option = "viridis")) )

      rm(loco)
      rm(htmxord)
      gc()

      subplot(hm_env$htpt1, margin = 0.01) %>%
        subplot(hm_env$htpt2,nrows=2, heights = c((1-input$hmgnprop),input$hmgnprop),shareX = TRUE,margin=0.01)
      }
    
    })
    
  }) # ,height=hm_heign2, width=hm_widz2)   MARCHE PAS ??? Try to find a way to resize or at least change pt1/pt2 proportions ?
  
  gc()
  ## How much memory is in use ?
  # output$foo5 <- renderText({
  #   req(input$click2) #to prevent print at first lauch
  #   isolate({
  #     memus<-mem_used() 
  #     print(paste0("memusage: ",memus))
  #   })
  # })
  output$topmodout <- renderText({
    print(paste0(gsub("corr_z_","",hm_env$mod_clu[as.numeric(input$topmod)]),collapse = ","))
    
  })
  ## Chi test to get genes differential between clusters [Not working anymore? Adapt to the matrix format ?]
  output$chirep<- renderText({
    req(input$click_run_chi) #to prevent print at first lauch
    isolate({
      clusterschi<-unlist(strsplit(input$clusterschi,split=","))                                                                          
      clusterschi<-gsub(" ","",clusterschi) 
      cel_in_chi<-c()
      for(i in as.numeric(clusterschi)){
        cel_in_chi<-append(cel_in_chi,unlist(hm_env$mcs_clu[i])) }
      cel_in_chi
      counts=dgcmtxcounts[,cel_in_chi]
      
      gene_mask=apply(counts,1,max)>3
      counts=counts[gene_mask,,drop=F]
      metacell_tot=colSums(counts)
      testarray<-array(c(counts,matrix(metacell_tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts))
      
      testchi<-suppressWarnings({apply(counts,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}) })
      testchi<-t(testchi)
      testchi=testchi[!is.na(testchi[,1]),]
      
      res=cbind(testchi,adjp=p.adjust(testchi[,1],method="BH"))
      mask=rownames(res)%in%res[,3]<0.01&  rowSums(counts)>=10
      
      a=res[mask,]
      genes_to_show=head(rownames(a)[order(a[,2],decreasing=T)],30)
      
      genes_ts_good=paste(genes_to_show,collapse=",")
      rm(counts)
      rm(testarray)
      print(paste0(genes_ts_good))
    })
  })
  
  output$mod_clust_bar <- renderPlot({
    print("plotting")
    df_zmod<-as.data.frame(hm_zmod)
    df_zmod$cluster<-hm_env$out$Cluster[match(rownames(df_zmod),hm_env$out$Metacell )]
    df_zmod$cluster<-factor(df_zmod$cluster, levels = paste0("cluster", 1:length(hm_env$metacells_kclust)))
    # 
    print(ggplot(data = df_zmod, aes(x = cluster, y = !!sym(paste0("My_Mod_",input$var2,"scaled")) )) +
            geom_boxplot(width=0.8, fill="white",size=0.35,outlier.alpha = 0.5)+
            xlab("Cluster")+theme_bw()+
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.85,size = rel(1)),
                  axis.title.x =element_text(size=rel(1),face="bold"),
                  axis.text.y = element_text(size = rel(1.5)),
                  axis.title.y =element_text(size=rel(1),face="bold"),
                  plot.title = element_text(size = rel(1.5), face = "bold"),
                  plot.subtitle = element_text(size = rel(1)),
                  panel.border = element_blank(),
                  legend.position = "none",
                  axis.line = element_line(colour = "black",size = 0.6)) + geom_jitter(height = 0,size=0.3, width = 0.2,alpha=0.3))
    
  })
  
  # Show all metacells selected in tab Mc annotation
  output$out_mcs = renderText({
    ht_mc_ord<-rownames(hm_zmod)[as.vector(unlist(column_order(ht))) ]  
    
    if(!input$ann_clu=="" ){
      ann_clu_show<-unlist(strsplit(input$ann_clu,split=","))
      ann_clu_show<-gsub(" ","",ann_clu_show)
      #Get the clusters from the heatmap (as cluster 1 is not named "1", I use numeric to get the first one, whatever its name)
      # And affect it to rownames of  hm_zmod as metacell "42" in heatmap is named "mc41" in the end
      mc_part1<-rownames(hm_zmod)[as.vector(unlist(  column_order(ht)[as.numeric(ann_clu_show)]   )) ]
    }else{mc_part1=""}
    
    
    if(!input$ann_add_mc1=="" && !input$ann_add_mc2=="" ){
      mc_part2<-c( ht_mc_ord[c(which(ht_mc_ord==paste0("mc",input$ann_add_mc1)):which(ht_mc_ord==paste0("mc",input$ann_add_mc2)))])
    }else{mc_part2=""}
    if(!input$ann_rem_mc1=="" && !input$ann_rem_mc2=="" ){
      mc_part3<-c( ht_mc_ord[c(which(ht_mc_ord==paste0("mc",input$ann_rem_mc1)):which(ht_mc_ord==paste0("mc",input$ann_rem_mc2)))])
    }else{mc_part3=""}
    mc_selected<-setdiff( append(mc_part1,mc_part2), mc_part3 )
    hm_env$current_mc_selected<-mc_selected
    print(mc_selected)
  })
  
  
  env_ann_df <- reactiveValues(df = ann_df)
  observeEvent(input$btn_save_name, {
    
    annotation <- input$grp_name
    
    if (!is.null(hm_env$current_mc_selected) && !is.null(input$grp_name)) {
      # Update the annotation in the dataframe
      env_ann_df$df[env_ann_df$df$metacell %in% hm_env$current_mc_selected, "annotation"] <- input$grp_name
    }
  })
  
  # Display cell proportions from given names
  output$out_prop_cellannot<- renderTable({
    table(env_ann_df$df$annotation)
  })
  # Save the annotations in a csv, give file name and choose where to save
  output$save_csv <- downloadHandler(
    filename = function() {
      paste0(input$file_name, ".csv")
    },
    content = function(file) {
      write.csv(env_ann_df$df, file, quote = F,col.names=T,row.names=F,sep = ",")
    })
  loadednames <- reactiveVal(NULL)
  
  observeEvent(input$annot_csv, {
    if (!is.null(input$annot_csv)) {
      loadednames(read.csv(input$annot_csv$datapath))
    }
  })
  
  observeEvent(input$btn_load_name, {
    env_ann_df$df <- loadednames()
  })
  output$reactive_txt <- renderUI({
    if (!input$name_fromlist == 'Gene (log CPM) or Gene signature') return(NULL) else {
      textInput("Not_on_list", "Write the gene or create a gene signature:")
    }
  })
  
  # Add or not the legend to the Bxplt
  legend_state <- reactiveVal(TRUE)
  observeEvent(input$toggle_legend, {
    legend_state(!legend_state())
    new_label <- ifelse(legend_state(), "Legend off", "Legend on")
    updateActionButton(session, "toggle_legend", label = new_label)
  })
  
  output$mod_group_bar <- renderPlot({
    print("plotting")
    # df_zmod<-as.data.frame(hm_zmod)
    #df_zmod$annot<-env_ann_df$df$annotation[match(rownames(df_zmod),env_ann_df$df$metacell )]
    df_zmod2$Annotation<-env_ann_df$df$annotation[match(rownames(df_zmod2),env_ann_df$df$metacell )]
    df_zmod2$currentscore<-NA
    if (input$name_fromlist == 'Gene (log CPM) or Gene signature'){
      var_select <- input$Not_on_list 
      ilist<-gene_input2vect(var_select,dgcmtxcounts)
      if(length(ilist)==1){
        eval(parse(text=paste0("df_zmod2$currentscore[match(rownames(df_zmod2),colnames(dgcmtxcounts))]<-preloco[ilist,,drop=F]")))
        
      }else if(length(ilist)>=2){
        
        values<- scale( Matrix::colSums(preloco_nc[ilist,,drop=F])/allsum )[,]  # Need [,] to drop the attributes of scale(x) so it works for after
        eval(parse(text=paste0("df_zmod2$currentscore[match(rownames(df_zmod2),names(values))]<-values")))
      }
      
      plot_w4<-ggplot(data = df_zmod2, aes(x = Annotation, y = currentscore, fill=Annotation, color=Annotation )) 
      
    }else {
         var_select <- input$name_fromlist
         plot_w4<-ggplot(data = df_zmod2, aes(x = Annotation, y = !!sym(paste0(var_select)),fill=Annotation, color=Annotation )) 
    }
    
    # Personalized palette 
    if(input$palette_1 == ""){
      palette_en_cours<-rep("black",50)
    }else{eval(parse(text=paste0("palette_en_cours<-",input$palette_1))) }
    # Personalized order
    if(input$bx_order == ""){
      bx_ord<-NULL
    }else{eval(parse(text=paste0("bx_ord<-",input$bx_order))) }
    # Show legend
    if (input$toggle_legend==T) {
      lgd_status="right"
    }else{lgd_status="none"}
    
    print(plot_w4+
            geom_boxplot(width=0.8, size=0.95*input$bx_siz,outlier.alpha = 0.5)+ #fill="white",
            xlab("Cluster")+theme_bw()+
            scale_fill_manual(values = alpha(palette_en_cours, .4)) +
            scale_color_manual(values = palette_en_cours) +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.85,size = rel(1)),
                  axis.title.x =element_text(size=rel(1),face="bold"),
                  #axis.text.y = element_text(size = rel(1.5)),
                  axis.text.y = element_text(size = rel(3)),
                  #axis.title.y =element_text(size=rel(1),face="bold"),
                  axis.title.y=element_blank(),
                  plot.title = element_text(size = rel(1.5), face = "bold"),
                  plot.subtitle = element_text(size = rel(1)),
                  panel.border = element_blank(),
                  legend.key.size = unit(1*input$bx_siz,'cm'),
                  legend.position = lgd_status,
                  legend.text = element_text(size=rel(1)*input$bx_siz),
                  axis.line = element_line(colour = "black",size = 0.6)) +
            geom_jitter(height = 0,size=0.5*input$bx_siz, width = 0.2,alpha=0.5) +
            scale_x_discrete(limits=bx_ord)
    )
    
  },height=bxhei, width=bxwid)
}


shinyApp(ui = ui, server = server)



#c('Fibro Os'='#C966EB','Fibro Bs'='#864036','Fibro Ds'='#D29FE3','Fibro Es'='#F87774','Fibro Ks'='#A3DDEB','Fibro Is'='#5C518A','bad'='#BDBDBD')
#c('Mono1'='#874037','Mono3'='#f77774','Mono4'='#ffa3a1','Mono5'='#B4DEC6','Macro6'='#c967eb','Macro7'='#5fc2ed','Macro8'='#4F88B9','Macro9'='#5C538B')
#c('Mono1','Mono3','Mono4','Mono5','Macro6','Macro7','Macro8','Macro9')

#c('Activated'='#8c1604','CD8s_MAIT_like'='#74BBCD','CM-like-1'='#e88935','CM-like-2'='#f0a869','CXCR6_TRM-like'='#cf9e9d','DP8alpha'='#e6a8da','Effector_CD8s'='#224f80','Effector_CD8s_GZMK_low'='#4F88B9','Effector_Tregs'='#7d4796','Tregs'='#be8ed4','ILC3-like'='#6e9160','NK-like'='#8ec977','gd-like_GZMB-'='#fbfca2','gd-like_GZMB+'='#c4de7e','Naive/CM'='#f0dac7','TBD_2'='#BCBCBB','TBD_4'='#999999','Tph-like'='#60918c','TRM-like-1'='#fc9795','TRM-like-2'='#f77774')
#c('Activated','CD8s_MAIT_like','CM-like-1','CM-like-2','CXCR6_TRM-like','DP8alpha','Effector_CD8s','Effector_CD8s_GZMK_low','Effector_Tregs','Tregs','ILC3-like','NK-like','gd-like_GZMB-','gd-like_GZMB+','Naive/CM','TBD_2','TBD_4','Tph-like','TRM-like-1','TRM-like-2')


#c('GC_like'='#ffd9d0','Atypical_Mem'='#e74c3c','Naive'='#ffad34','Mem_1'='#4575BB','Mem_2_NFkB'='#004b8d','Mem_3'='#F09fB6')
#c('Naive','Mem_1','Mem_2_NFkB','Mem_3','Atypical_Mem','GC_like')
# 31,30,29,27,28,18,24,23,25,20,21,22,26,16,17,15,32,39,1,33,37,36,38,8,9,12,11,10,7,3,2,13,35,4,6,14,19,5,40,(34)
