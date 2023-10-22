library(reshape2)
library(stringi)
library(pbapply)
library(ggplot2)
library(sf)
library(raster)
library(parallel)
library(ggpubr)
library(plotrix)
library(viridis)
library(stringr)

##### Reproductive system ----
sisrep= read.csv("datasets/sistreprodutivo.csv")
sisrep

sisrep_noNA=na.omit(sisrep)
x2_data<-table(sisrep_noNA$origen, sisrep_noNA$sistrep)
chisq.test(x2_data)

library(reshape2)

melt(t(x2_data))

data1= melt(t(x2_data))

library(ggplot2)
library(viridis)
#devtools::install_github("BlakeRMills/MoMAColors")
library(MoMAColors)

sistreplot= ggplot(data1, aes(y=value, x= Var2, fill= Var1))+
  geom_bar(position= "fill", stat= "identity")+
  coord_flip()+
  scale_fill_manual(labels= c("Dioecious", "Self-compatible", "Self-incompatible"),values=moma.colors("VanGogh", 3), name= "Reproductive\nsystem")+
  scale_x_discrete(name= "Plant origin", labels= c("Exotic", "Native"))+
  scale_y_continuous(name= "Proportion of crop species (%)")+
  theme_bw()+
  theme(legend.position = "right",legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=14), legend.text = element_text(family ="sans", size=12), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));sistreplot

##### Spatial distribution of crops -----

poli= read.csv("datasets/category.csv",row.names = 1)
poli= poli[-nrow(poli),]
poli$Latin.name=word(poli$Latin.name, 1, 2)
data2= data.frame(table(poli$Latin.name, poli$group))
rich= aggregate(data2$Freq, list(data2$Var1, data2$Var2), sum)
rich$origin= poli[match(rich$Group.1, poli$Latin.name),4]
rich
colnames(rich)= c("spp", "group", "rich", "origin")

aov1= aov(rich~origin*group, data= rich)
summary(aov1)

as.data.frame.summary.aovlist <- function(x) {
  if(length(x) == 1) {
    as.data.frame(x[[1]])
  } else {
    lapply(unlist(x, FALSE), as.data.frame)
  }
}
write.csv(as.data.frame.summary.aovlist(summary(aov1)),"results/pol_richness_AOV.csv")

TukeyHSD(aov(rich~group*origin,data=rich))
write.csv(TukeyHSD(aov(rich~group*origin,data=rich))$group,"results/pol_richness_Tukey_polgroup.csv")
write.csv(TukeyHSD(aov(rich~group*origin,data=rich))$origin,"results/pol_richness_Tukey_croporigin.csv")
write.csv(TukeyHSD(aov(rich~group*origin,data=rich))$`group:origin`,"results/pol_richness_Tukey_interaction.csv")

library(ggplot2)
library(viridis)
library(plotrix)

aggregate(rich$rich,list(rich$group,rich$origin),FUN = function(x) c(mean = mean(x), se = std.error(x)))->Meanz
colnames(Meanz)[1:2]=c("group","origin")

Meanz$lower<-Meanz$x[,1]-(Meanz$x[,2]*2)
Meanz$upper<-Meanz$x[,1]+(Meanz$x[,2]*2)

poli_group= ggplot(Meanz, aes(x=factor(origin,levels=c("Nativa","Exotica")), y= x[,1], color= origin))+
  facet_wrap(.~group,scales = "free")+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,linewidth=.6) +
  geom_point(size=2)+
  scale_x_discrete(labels=c("Native","Exotic"),name="Crop origin")+
  scale_color_manual(values=c("#EB5353","#36AE7C"))+
  scale_y_continuous(name= "Richness")+
  theme_bw()+
  theme(strip.text = element_text(size=10,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); poli_group

tiff("raw_figures/Figure S1.tiff", units= "in", width= 4.5, height= 3.5, res= 900)
plot(poli_group)
dev.off()

## Pollinator Webs ----
library(ggsankey)
AllData=table(poli$Latin.name, poli$pollinator.spp)
othergroupsAll<-{}
for(x in colnames(AllData)){
  othergroupsAll[x]=unique(poli[grep(x,poli$pollinator.spp),]$group)
}

groupsplanta_All= rownames(AllData)[order(rownames(AllData))]
data.frame(Plant=groupsplanta_All,Origin=NA)->ExotAndNat
for(x in 1:nrow(ExotAndNat)){
  unique(poli[grep(ExotAndNat[x,1],poli$Latin.name),]$Brasil)->ExotAndNat[x,2]
}
PlantExo<-ExotAndNat[ExotAndNat$Origin=="Exotica",1]
PlantNat<-ExotAndNat[ExotAndNat$Origin=="Nativa",1]
groupsbees_All= colnames(AllData)[othergroupsAll=="Bees"][order(colnames(AllData)[othergroupsAll=="Bees"])]
groupsIns_All= colnames(AllData)[othergroupsAll=="Other insects"][order(colnames(AllData)[othergroupsAll=="Other insects"])]
groupsVert_All= colnames(AllData)[othergroupsAll=="Vertebrates"][order(colnames(AllData)[othergroupsAll=="Vertebrates"])]

reshape2::melt(AllData)->IntAll

dfAll <- IntAll[IntAll$value==1,] %>%
  make_long(Var2,Var1)
dfAll$col<-NA
for(x in 1:nrow(dfAll)){
  
  FoundPltExo<-grep(paste0("^",dfAll$node[x],"$"),PlantExo)
  FoundPltNat<-grep(paste0("^",dfAll$node[x],"$"),PlantNat)
  FoundBee<-grep(paste0("^",dfAll$node[x],"$"),groupsbees_All)
  FoundIns<-grep(paste0("^",dfAll$node[x],"$"),groupsIns_All)
  FoundVert<-grep(paste0("^",dfAll$node[x],"$"),groupsVert_All)
  
  if(!is.na(FoundPltExo[1])){
    dfAll$col[x]<-"Exotic plant"  
  }
  if(!is.na(FoundPltNat[1])){
    dfAll$col[x]<-"Native plant"  
  }
  if(!is.na(FoundBee[1])){
    dfAll$col[x]<-"Bees"  
  }
  if(!is.na(FoundIns[1])){
    dfAll$col[x]<-"Other insects"  
  }
  if(!is.na(FoundVert[1])){
    dfAll$col[x]<-"Vertebrates"  
  }
}

WebAll<-ggplot(dfAll, aes(x = x, 
                          next_x = next_x, 
                          node = factor(node,levels=c(groupsbees_All,groupsIns_All,PlantExo,PlantNat,groupsVert_All)), 
                          fill=factor(col,levels=c("Bees","Exotic plant","Native plant","Other insects","Vertebrates")),
                          next_node = next_node)) +
  geom_sankey(flow.alpha = 0.75, node.color = 1,linewidth=0.01) +
  coord_flip()+
  scale_fill_manual(name="Group",labels=c("Bees","Exotic plants", "Native plants", "Other insects", "Vertebrates"),values=c("#F9D923","#EB5353","#36AE7C","orange","#187498"))+
  theme_void()+
  theme(legend.text.align = 0,legend.title.align=0,legend.position = "right",legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10));WebAll

ggplot_build(WebAll)->BuiltAll

SelNames<-(BuiltAll$data[[2]]$ymax-BuiltAll$data[[2]]$ymin)>quantile(BuiltAll$data[[2]]$ymax-BuiltAll$data[[2]]$ymin)[4]

data.frame(node=as.character(BuiltAll$data[[2]]$node),group=NA)->UniquesAll
for(x in 1:nrow(UniquesAll)){
  FoundExo<-grep(paste0("^",UniquesAll$node[x],"$"),PlantExo)
  FoundNat<-grep(paste0("^",UniquesAll$node[x],"$"),PlantNat)
  FoundBee<-grep(paste0("^",UniquesAll$node[x],"$"),groupsbees_All)
  FoundIns<-grep(paste0("^",UniquesAll$node[x],"$"),groupsIns_All)
  FoundVert<-grep(paste0("^",UniquesAll$node[x],"$"),groupsVert_All)
  
  if(!is.na(FoundExo[1])){
    UniquesAll$group[x]<-"Exotic plant"  
  }
  if(!is.na(FoundNat[1])){
    UniquesAll$group[x]<-"Native plant"  
  }
  if(!is.na(FoundBee[1])){
    UniquesAll$group[x]<-"Bees"  
  }
  if(!is.na(FoundIns[1])){
    UniquesAll$group[x]<-"Other insects"  
  }
  if(!is.na(FoundVert[1])){
    UniquesAll$group[x]<-"Vertebrates"  
  }
}
PltAll<-UniquesAll
PltAll$node<-paste0("         ",substr(PltAll$node, 1, 1), ". ", word(PltAll$node, 2, 2))
PltAll[!SelNames,]<-" "
PltAll[!grepl("plant",PltAll$group),]$node<-""

WebAll<-WebAll+geom_sankey_text(label=PltAll$node, size = 2, color = "black",angle=90,hjust = 0,fontface="italic")

InsAll<-UniquesAll
InsAll$node<-paste0(substr(InsAll$node, 1, 1), ". ", word(InsAll$node, 2, 2),"         ")
InsAll[!SelNames,]<-" "
InsAll[grepl("plant",PltAll$group),][,1]<-" "

WebAll<-WebAll+geom_sankey_text(label=InsAll$node, size = 2, color = "black",angle=90,hjust = 1,fontface="italic"); WebAll

tiff("raw_figures/Figure 1.tiff", units= "in", width= 14, height= 7, res= 900)
WebAll
dev.off()

##### Maps ----

read.csv("raw_datasets/PAM - 2021/area_colhida.csv",sep=";")->area_col
read.csv("raw_datasets/PAM - 2021/area_plantada.csv",sep=";")->area_plt
read.csv("raw_datasets/PAM - 2021/quant.csv",sep=";")->quant
read.csv("raw_datasets/PAM - 2021/rend.csv",sep=";")->rend
read.csv("raw_datasets/PAM - 2021/valor.csv",sep=";")->valor

for(x in 2:ncol(area_col)){
  gsub("-",0,area_col[,x])->area_col[,x]
  gsub("-",0,area_plt[,x])->area_plt[,x]
  gsub("-",0,quant[,x])->quant[,x]
  gsub("-",0,rend[,x])->rend[,x]
  gsub("-",0,valor[,x])->valor[,x]
}

for(x in 2:ncol(area_col)){
  as.numeric(area_col[,x])->area_col[,x]
  as.numeric(area_plt[,x])->area_plt[,x]
  as.numeric(quant[,x])->quant[,x]
  as.numeric(rend[,x])->rend[,x]
  as.numeric(valor[,x])->valor[,x]
}

melt(area_col)->AllMuns
melt(area_plt)[,3]->AllMuns$AreaPlt
melt(quant)[,3]->AllMuns$Quant
melt(rend)[,3]->AllMuns$Rend
melt(valor)[,3]->AllMuns$Valor

AllMuns$Valor*1000->AllMuns$Valor

gsub("[\\(\\)]", "", regmatches(AllMuns$Município, gregexpr("\\(.*?\\)", AllMuns$Município)))->AllMuns$UF
colnames(AllMuns)[2:3]=c("Crop","AreaColh")
stri_trans_general(str = AllMuns$Crop,id = "Latin-ASCII")->AllMuns$Crop

AllMuns[!grepl("Borracha|Linho|Algodao|Tabaco|Fumo|Sisal|Rami|Tungue|Cana.para.forragem|Cafe..em.grao..Total",AllMuns$Crop),]->AllMuns

read.csv("datasets/final_dataset_crops.csv",sep=";",encoding = "latin1")->pol_depend
pol_depend$pol_depend_number<-pol_depend$pol_dep
gsub("no increase",0,pol_depend$pol_depend_number)->pol_depend$pol_depend_number
gsub("little",0.25,pol_depend$pol_depend_number)->pol_depend$pol_depend_number
gsub("modest",0.5,pol_depend$pol_depend_number)->pol_depend$pol_depend_number
gsub("great",0.75,pol_depend$pol_depend_number)->pol_depend$pol_depend_number
gsub("essential",1,pol_depend$pol_depend_number)->pol_depend$pol_depend_number
as.numeric(pol_depend$pol_depend_number)->pol_depend$pol_depend_number

read.csv("datasets/crops_info.csv",sep=";",encoding = "latin1")->crops_info
crops_info$origin<-pol_depend[match(crops_info$spp,pol_depend$spp),]$origin
crops_info$pol_dep<-pol_depend[match(crops_info$spp,pol_depend$spp),]$pol_depend_number

AllMuns$crop_origin<-crops_info[match(AllMuns$Crop,crops_info$shape_name),]$origin

EVP<-pbapply::pblapply(1:nrow(AllMuns), function(x) {
  (AllMuns$Valor[x]*crops_info[grep(paste0("^",AllMuns$Crop[x],"$"),crops_info$shape_name),]$pol_dep)
})
unlist(EVP)->AllMuns$EVP

AllMuns_NoZero=AllMuns[!AllMuns$Valor==0,]
AllMuns_NoZero$NetEVP<-vegan::decostand(abs((AllMuns_NoZero$Valor-AllMuns_NoZero$EVP)-max(AllMuns_NoZero$Valor-AllMuns_NoZero$EVP,na.rm = T)),method="range",na.rm=T)

(sum(AllMuns_NoZero$Valor,na.rm = T)*0.1855) #Total earnings of agriculture in Brazil in 2021. Total of $ 131,352,000,000
(sum(AllMuns_NoZero$EVP,na.rm = T)*0.1855) #Economic value of Pollinator in US dollars based on the average exchange rate from 2021. Total of $41,458,757,637
aggregate(AllMuns_NoZero$EVP,list(AllMuns_NoZero$Crop),sum)->Earn_Crops
rbind(Earn_Crops,data.frame(Group.1="Café total",x=sum(Earn_Crops[grepl("Cafe",Earn_Crops$Group.1),"x"])))->Earn_Crops
Earn_Crops[!grepl("Cafe",Earn_Crops$Group.1),]->Earn_Crops
Earn_Crops[order(Earn_Crops$x,decreasing = T),][1:5,]->Top5
Top5$x<-Top5$x*0.1855
Top5 # Economic value of the five crops that contributed more to economic earnings. Soybean = $31,697,095,465. Cotton = $4,917,113,240. Total coffee = $3,236,660,670. Orange = $1,162,600,938. Açai = $738,132,279
sum(Top5$x)/(sum(AllMuns_NoZero$EVP,na.rm = T)*0.1855)

aggregate(AllMuns$EVP,list(AllMuns$Município),sum,na.rm=T)->EVPMun
aggregate(AllMuns[!grepl("soja",AllMuns$Crop,ignore.case=T),]$EVP,list(AllMuns[!grepl("soja",AllMuns$Crop,ignore.case=T),]$Município),sum,na.rm=T)->EVPMun_NoSoy

aggregate(AllMuns$Valor,list(AllMuns$Município),sum,na.rm=T)->ValorAll
aggregate(AllMuns[!grepl("soja",AllMuns$Crop,ignore.case=T),]$Valor,list(AllMuns[!grepl("soja",AllMuns$Crop,ignore.case=T),]$Município),sum,na.rm=T)->Valor_NoSoy

EVPMun$x*1/ValorAll$x->EVPMun$NetEVP
EVPMun_NoSoy$x*1/Valor_NoSoy$x->EVPMun_NoSoy$NetEVP

AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Exo",]->AllMuns_Exo
(sum(AllMuns_Exo$EVP,na.rm=T)*0.1855) #Economic value of pollination in US dollars for Exotic crops based on the average exchange rate from 2021. Total of $38,571,729,036
aggregate(AllMuns_Exo$EVP,list(AllMuns_Exo$Município),sum,na.rm=T)->EVPMun_Exo
aggregate(AllMuns_Exo[!grepl("soja",AllMuns_Exo$Crop,ignore.case=T),]$EVP,list(AllMuns_Exo[!grepl("soja",AllMuns_Exo$Crop,ignore.case=T),]$Município),sum,na.rm=T)->EVPMun_Exo_NoSoy
aggregate(AllMuns_Exo$Valor,list(AllMuns_Exo$Município),sum,na.rm=T)->ValorAll_Exo
aggregate(AllMuns_Exo[!grepl("soja",AllMuns_Exo$Crop,ignore.case=T),]$Valor,list(AllMuns_Exo[!grepl("soja",AllMuns_Exo$Crop,ignore.case=T),]$Município),sum,na.rm=T)->ValorAll_Exo_NoSoy
EVPMun_Exo$x*1/ValorAll_Exo$x->EVPMun_Exo$NetEVP
EVPMun_Exo_NoSoy$x*1/ValorAll_Exo_NoSoy$x->EVPMun_Exo_NoSoy$NetEVP

AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Nat",]->AllMuns_Nat
(sum(AllMuns_Nat$EVP,na.rm = T)*0.1855) #Economic value of pollination in US dollars for Native crops based on the average exchange rate from 2021. Total of $2,887,028,601
aggregate(AllMuns_Nat$EVP,list(AllMuns_Nat$Município),sum,na.rm=T)->EVPMun_Nat
aggregate(AllMuns_Nat$Valor,list(AllMuns_Nat$Município),sum,na.rm=T)->ValorAll_Nat
EVPMun_Nat$x*1/ValorAll_Nat$x->EVPMun_Nat$NetEVP

na.omit(AllMuns_NoZero)->AllMuns_NoZero
data.frame(rowSums(table(AllMuns_NoZero$Município,AllMuns_NoZero$Crop)))->Richness
colnames(Richness)="S"
data.frame(rowSums(table(AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Nat",]$Município,AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Nat",]$Crop)))->Richness_Nat
data.frame(rowSums(table(AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Exo",]$Município,AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Exo",]$Crop)))->Richness_Exo
data.frame(rowSums(table(AllMuns_NoZero$Município,AllMuns_NoZero$Crop)))->Richness
colnames(Richness)="S"
colnames(Richness_Exo)="S"
colnames(Richness_Nat)="S"

aggregate(AllMuns_NoZero$AreaPlt,list(AllMuns_NoZero$Município),sum)->Area_All
aggregate(AllMuns_NoZero[!grepl("Soja",AllMuns_NoZero$Crop),]$AreaPlt,list(AllMuns_NoZero[!grepl("Soja",AllMuns_NoZero$Crop),]$Município,AllMuns_NoZero[!grepl("Soja",AllMuns_NoZero$Crop),]$crop_origin),sum)->AreaOriginNoSoy
AreaOriginNoSoy[AreaOriginNoSoy$Group.2=="Exo",]->AreaOriginNoSoy
AreaOriginNoSoy$Group.2<-"Exo_NoSoy"
aggregate(AllMuns_NoZero$AreaPlt,list(AllMuns_NoZero$Município,AllMuns_NoZero$crop_origin),sum)->Area_Origin
rbind(AreaOriginNoSoy,Area_Origin)->Area_Origin

st_read("datasets/shapefiles/BR_Municipios_2022.shp")->BRMun
BRMun$EVP<-NA
BRMun$ValorTot<-NA
BRMun$NetEVP<-NA
BRMun$S<-NA
BRMun$S_Nat<-NA
BRMun$S_Exo<-NA
BRMun$NetEVP_Nat<-NA
BRMun$NetEVP_Exo<-NA
BRMun$NetEVP_Exo_NoSoy<-NA
BRMun$NetEVP_All_NoSoy<-NA
BRMun$Area_All<-NA
BRMun$Area_Nat<-NA
BRMun$Area_Exo<-NA
BRMun$Area_ExoNoSoy<-NA


BRMun[match(Area_All$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$Area_All<-Area_All$x
BRMun[match(Area_Origin[Area_Origin$Group.2=="Nat",]$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$Area_Nat<-Area_Origin[Area_Origin$Group.2=="Nat",]$x
BRMun[match(Area_Origin[Area_Origin$Group.2=="Exo",]$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$Area_Exo<-Area_Origin[Area_Origin$Group.2=="Exo",]$x
BRMun[match(Area_Origin[Area_Origin$Group.2=="Exo_NoSoy",]$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$Area_ExoNoSoy<-Area_Origin[Area_Origin$Group.2=="Exo_NoSoy",]$x
BRMun[match(EVPMun$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$EVP<-EVPMun$x
BRMun[match(ValorAll$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$ValorTot<-ValorAll$x
BRMun[match(EVPMun$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$NetEVP<-EVPMun$NetEVP
BRMun[match(rownames(Richness),paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$S<-Richness$S
BRMun[match(rownames(Richness_Nat),paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$S_Nat<-Richness_Nat$S
BRMun[match(rownames(Richness_Exo),paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$S_Exo<-Richness_Exo$S
BRMun[match(EVPMun_Nat$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$NetEVP_Nat<-EVPMun_Nat$NetEVP
BRMun[match(EVPMun_Exo$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$NetEVP_Exo<-EVPMun_Exo$NetEVP
BRMun[match(EVPMun_NoSoy$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$NetEVP_All_NoSoy<-EVPMun_NoSoy$NetEVP
BRMun[match(EVPMun_Exo_NoSoy$Group.1,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$NetEVP_Exo_NoSoy<-EVPMun_Exo_NoSoy$NetEVP

st_read("datasets/shapefiles/dashboard_biomes-static-layer.shp",options=c("encoding=latin1"))->Biomes
st_write(BRMun,"FinalMunInfo_Final.shp",append=FALSE)
#st_read("FinalMunInfo_Final.shp")->BRMun
#colnames(BRMun)=c("CD_MUN","NM_MUN","SIGLA_UF","AREA_KM","EVP","ValorTot","NetEVP","NetEVP_Nat","NetEVP_Exo","NatVeg30Km","NatVegLoc","PolRiskAll","PolRiskExo","PolRiskNat","S","S_Nat","S_Exo","NetEVP_Exo_NoSoy","NetEVP_All_NoSoy","PolRiskAll_NoSoy","PolRiskExo_NoSoy","Area_All","Area_Nat","Area_Exo","Area_ExoNoSoy","geometry")

##### Maps for richness ----
S_All<-ggplot() +
  geom_sf(data = BRMun,aes(color=S,fill=S)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "#362FD9",name="Crop\nrichness")+
  scale_fill_gradient(low = "white",high = "#362FD9",name="Crop\nrichness")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

S_Nat<-ggplot() +
  geom_sf(data = BRMun,aes(color=S_Nat,fill=S_Nat)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "#06FF00",name="Crop\nrichness")+
  scale_fill_gradient(low = "white",high = "#06FF00",name="Crop\nrichness")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

S_Exo<-ggplot() +
  geom_sf(data = BRMun,aes(color=S_Exo,fill=S_Exo)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "#FF1700",name="Crop\nrichness")+
  scale_fill_gradient(low = "white",high = "#FF1700",name="Crop\nrichness")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

tiff("raw_figures/Richness_MainMaps - Part I.tiff", units="in", width=10, height=5, res=900) #Save plot
ggarrange(S_Nat,S_Exo,nrow = 1,ncol = 2)
dev.off()

##### Maps of EVP ----

PolDepMap_All<-ggplot() +
  geom_sf(data = BRMun,aes(color=NetEVP,fill=NetEVP)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_fill_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolDepMap_All_NoSoy<-ggplot() +
  geom_sf(data = BRMun,aes(color=NetEVP_All_NoSoy,fill=NetEVP_All_NoSoy)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_fill_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolDepMap_Exo<-ggplot() +
  geom_sf(data = BRMun,aes(color=NetEVP_Exo,fill=NetEVP_Exo)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_fill_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolDepMap_Exo_NoSoy<-ggplot() +
  geom_sf(data = BRMun,aes(color=NetEVP_Exo_NoSoy,fill=NetEVP_Exo_NoSoy)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_fill_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolDepMap_Nat<-ggplot() +
  geom_sf(data = BRMun,aes(color=NetEVP_Nat,fill=NetEVP_Nat)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_fill_gradient(low = "white",high = "green",name="Pollinator\ndependency")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

tiff("raw_figures/PolDep_MainMaps - Part I.tiff", units="in", width=12, height=5, res=900) #Save plot
ggarrange(PolDepMap_Nat,PolDepMap_Exo,PolDepMap_Exo_NoSoy,common.legend = TRUE,legend = "right",nrow = 1,ncol = 3)
dev.off()

colnames(EVPMun_Nat)=c("Mun","EVP","NetEVP")
colnames(EVPMun_Exo)=c("Mun","EVP","NetEVP")
rbind(data.frame(EVPMun_Nat,origin="Native"),data.frame(EVPMun_Exo,origin="Exotic"))->NetEVPAll

###### Calculating vegetation -----

raster("datasets/rasters/brasil_coverage_2021.tif")->Brazil_cov
BRMun$NatVeg30Km<-NA
BRMun$NatVegLoc<-NA

#library(parallel)
#cl <- makeCluster(rep("localhost",4), type="SOCK")
#clusterExport(cl = cl,list("BRMun","Brazil_cov"))
#clusterEvalQ(cl,library(sf))
#clusterEvalQ(cl,library(raster))

NatVeg_All<-pbapply::pblapply(1:nrow(BRMun), function(x) {  
  try(st_buffer(BRMun[x,],30000))->Buf30km
  if(is(Buf30km)[1]=="try-error"){
    return(list(NA,NA))
    next
  }
  r2 <- try(raster::crop(Brazil_cov, raster::extent(Buf30km))) ### crop to the extent
  if(is(r2)[1]=="try-error"){
    return(list(NA,NA))
    next
  }
  maskalt = raster::mask(r2, Buf30km)
  AllValues<-raster::getValues(maskalt)
  rm(Buf30km,r2,maskalt)
  gc()
  data.frame(table(na.omit(AllValues)))->Freq ## 1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 50, 13 are all native vegetation
  Freq[!Freq[,1]==0,]->Freq
  NatVeg<-Freq[grepl(paste0("^",c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 50, 13),"$",collapse = "|"),Freq[,1]),]
  OtherVeg<-Freq[!grepl(paste0("^",c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 50, 13),"$",collapse = "|"),Freq[,1]),]
  PercNat30km<-sum(NatVeg[,2])*1/sum(NatVeg[,2],OtherVeg[,2])
  
  r2 <- raster::crop(Brazil_cov, raster::extent(BRMun[x,])) ### crop to the extent
  maskalt = raster::mask(r2, BRMun[x,])
  rm(r2)
  gc()
  AllValues<-raster::getValues(maskalt)
  rm(maskalt)
  gc()
  data.frame(table(na.omit(AllValues)))->Freq ## 1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 50, 13 are all native vegetation
  Freq[!Freq[,1]==0,]->Freq
  NatVeg<-Freq[grepl(paste0("^",c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 50, 13),"$",collapse = "|"),Freq[,1]),]
  OtherVeg<-Freq[!grepl(paste0("^",c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 50, 13),"$",collapse = "|"),Freq[,1]),]
  PercNatLoc<-sum(NatVeg[,2])*1/sum(NatVeg[,2],OtherVeg[,2])
  
  return(list(PercNat30km,PercNatLoc))
})

#stopCluster(cl)
#save.image("D:/Luke Files/Ciencia/Papers/The one about plant-polinator network - Will et al (2023)/PolEco/.RData")

All_PercNat30km<-unlist(NatVeg_All)[seq(1,length(unlist(NatVeg_All)),2)]
All_PercNatLoc<-unlist(NatVeg_All)[seq(2,length(unlist(NatVeg_All)),2)]
BRMun$NatVeg30Km<-All_PercNat30km
BRMun$NatVegLoc<-All_PercNatLoc

BRMun$PolRiskAll<-NA
BRMun$PolRiskAll_NoSoy<-NA
BRMun$PolRiskExo<-NA
BRMun$PolRiskExo_NoSoy<-NA
BRMun$PolRiskNat<-NA

BRMun$PolRiskAll<-(BRMun$NetEVP*abs(BRMun$NatVeg30Km-1))#/2
BRMun$PolRiskAll_NoSoy<-(BRMun$NetEVP_All_NoSoy*abs(BRMun$NatVeg30Km-1))#/2
BRMun$PolRiskExo<-(BRMun$NetEVP_Exo*abs(BRMun$NatVeg30Km-1))#/2
BRMun$PolRiskExo_NoSoy<-(BRMun$NetEVP_Exo_NoSoy*abs(BRMun$NatVeg30Km-1))#/2
BRMun$PolRiskNat<-(BRMun$NetEVP_Nat*abs(BRMun$NatVeg30Km-1))#/2

st_write(BRMun,"datasets/shapefiles/FinalMunInfo.shp",append=FALSE)

##### Maps of Polination risk -----

PolRiskMap_All<-ggplot() +
  geom_sf(data = BRMun,aes(color=PolRiskAll,fill=PolRiskAll)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_fill_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolRiskMap_All_NoSoy<-ggplot() +
  geom_sf(data = BRMun,aes(color=PolRiskAll_NoSoy,fill=PolRiskAll_NoSoy)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_fill_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolRiskMap_Nat<-ggplot() +
  geom_sf(data = BRMun,aes(color=PolRiskNat,fill=PolRiskNat)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  scale_fill_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolRiskMap_Exo<-ggplot() +
  geom_sf(data = BRMun,aes(color=PolRiskExo,fill=PolRiskExo)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_fill_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

PolRiskMap_Exo_NoSoy<-ggplot() +
  geom_sf(data = BRMun,aes(color=PolRiskExo_NoSoy,fill=PolRiskExo_NoSoy)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_fill_gradient(low = "white",high = "red",name="Pollinator\nshortage risk")+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

tiff("raw_figures/PolRisk_MainMaps - Part II.tiff", units="in", width=12, height=5, res=900) #Save plot
ggarrange(PolRiskMap_Nat,PolRiskMap_Exo,PolRiskMap_Exo_NoSoy,common.legend = TRUE,legend = "right",nrow = 1,ncol = 3)
dev.off()

tiff("raw_figures/All_MainMaps.tiff", units="in", width=15, height=8, res=900) #Save plot
ggarrange(S_All,PolDepMap_All,PolDepMap_All_NoSoy,PolRiskMap_All,PolRiskMap_All_NoSoy,nrow = 2,ncol = 3)
dev.off()

##### For biomes ----
sf_use_s2(FALSE)
st_make_valid(BRMun)->BRMun_valid
st_make_valid(Biomes)->Biomes_valid
st_intersection(x = Biomes_valid[Biomes_valid$name=="Amazônia",],y = BRMun_valid)->Amz
st_intersection(x = Biomes_valid[Biomes_valid$name=="Caatinga",],y = BRMun_valid)->Caa
st_intersection(x = Biomes_valid[Biomes_valid$name=="Cerrado",],y = BRMun_valid)->Cer
st_intersection(x = Biomes_valid[Biomes_valid$name=="Mata Atlântica",],y = BRMun_valid)->Atl
st_intersection(x = Biomes_valid[Biomes_valid$name=="Pampa",],y = BRMun_valid)->Pmp
st_intersection(x = Biomes_valid[Biomes_valid$name=="Pantanal",],y = BRMun_valid)->Pnt

data.frame(Biome=c(Amz$name,Caa$name,Cer$name,Atl$name,Pmp$name,Pnt$name),
           Mun=c(Amz$NM_MUN,Caa$NM_MUN,Cer$NM_MUN,Atl$NM_MUN,Pmp$NM_MUN,Pnt$NM_MUN),
           State=c(Amz$SIGLA_U,Caa$SIGLA_U,Cer$SIGLA_U,Atl$SIGLA_U,Pmp$SIGLA_U,Pnt$SIGLA_U),
           NatVeg30km=c(Amz$NatVeg30Km,Caa$NatVeg30Km,Cer$NatVeg30Km,Atl$NatVeg30Km,Pmp$NatVeg30Km,Pnt$NatVeg30Km),
           PolDep_All=c(Amz$NetEVP,Caa$NetEVP,Cer$NetEVP,Atl$NetEVP,Pmp$NetEVP,Pnt$NetEVP),
           PolDep_All_NoSoy=c(Amz$NetEVP_All_NoSoy,Caa$NetEVP_All_NoSoy,Cer$NetEVP_All_NoSoy,Atl$NetEVP_All_NoSoy,Pmp$NetEVP_All_NoSoy,Pnt$NetEVP_All_NoSoy),
           PolDep_Nat=c(Amz$NetEVP_Nat,Caa$NetEVP_Nat,Cer$NetEVP_Nat,Atl$NetEVP_Nat,Pmp$NetEVP_Nat,Pnt$NetEVP_Nat),
           PolDep_Exo=c(Amz$NetEVP_Exo,Caa$NetEVP_Exo,Cer$NetEVP_Exo,Atl$NetEVP_Exo,Pmp$NetEVP_Exo,Pnt$NetEVP_Exo),
           PolDep_Exo_NoSoy=c(Amz$NetEVP_Exo_NoSoy,Caa$NetEVP_Exo_NoSoy,Cer$NetEVP_Exo_NoSoy,Atl$NetEVP_Exo_NoSoy,Pmp$NetEVP_Exo_NoSoy,Pnt$NetEVP_Exo_NoSoy),
           PolRisk_All=c(Amz$PolRiskAll,Caa$PolRiskAll,Cer$PolRiskAll,Atl$PolRiskAll,Pmp$PolRiskAll,Pnt$PolRiskAll),
           PolRisk_All_NoSoy=c(Amz$PolRiskAll_NoSoy,Caa$PolRiskAll_NoSoy,Cer$PolRiskAll_NoSoy,Atl$PolRiskAll_NoSoy,Pmp$PolRiskAll_NoSoy,Pnt$PolRiskAll_NoSoy),
           PolRisk_Nat=c(Amz$PolRiskNat,Caa$PolRiskNat,Cer$PolRiskNat,Atl$PolRiskNat,Pmp$PolRiskNat,Pnt$PolRiskNat),
           PolRisk_Exo=c(Amz$PolRiskExo,Caa$PolRiskExo,Cer$PolRiskExo,Atl$PolRiskExo,Pmp$PolRiskExo,Pnt$PolRiskExo),
           PolRisk_Exo_NoSoy=c(Amz$PolRiskExo_NoSoy,Caa$PolRiskExo_NoSoy,Cer$PolRiskExo_NoSoy,Atl$PolRiskExo_NoSoy,Pmp$PolRiskExo_NoSoy,Pnt$PolRiskExo_NoSoy),
           S_All=c(Amz$S,Caa$S,Cer$S,Atl$S,Pmp$S,Pnt$S),
           S_Nat=c(Amz$S_Nat,Caa$S_Nat,Cer$S_Nat,Atl$S_Nat,Pmp$S_Nat,Pnt$S_Nat),
           S_Exo=c(Amz$S_Exo,Caa$S_Exo,Cer$S_Exo,Atl$S_Exo,Pmp$S_Exo,Pnt$S_Exo),
           Area_All=c(Amz$Area_All,Caa$Area_All,Cer$Area_All,Atl$Area_All,Pmp$Area_All,Pnt$Area_All),
           Area_Nat=c(Amz$Area_Nat,Caa$Area_Nat,Cer$Area_Nat,Atl$Area_Nat,Pmp$Area_Nat,Pnt$Area_Nat),
           Area_Exo=c(Amz$Area_Exo,Caa$Area_Exo,Cer$Area_Exo,Atl$Area_Exo,Pmp$Area_Exo,Pnt$Area_Exo),
           Area_ExoNoSoy=c(Amz$Area_ExoNoSoy,Caa$Area_ExoNoSoy,Cer$Area_ExoNoSoy,Atl$Area_ExoNoSoy,Pmp$Area_ExoNoSoy,Pnt$Area_ExoNoSoy))->BiomesData

rm(Amz,Caa,Cer,Atl,Pmp,Pnt)
gc()
#write.csv(BiomesData,"BiomesData.csv")

##### Differences in Richness ----
melt(BiomesData)->BiomesTest
BiomesTest[!grepl("All",BiomesTest$variable),]->BiomesTest
BiomesTest[grepl("S_",BiomesTest$variable),]->S_Test
gsub("S_Nat","Native",S_Test$variable)->S_Test$variable
gsub("S_Exo","Exotic",S_Test$variable)->S_Test$variable
colnames(S_Test)[4:5]=c("Origin","Richness")
gsub("Amazônia","Amazon",S_Test$Biome)->S_Test$Biome
gsub("Pampa","Pampas",S_Test$Biome)->S_Test$Biome
gsub("Mata Atlântica","Atl. forest",S_Test$Biome)->S_Test$Biome

aov(Richness~Biome*Origin,data = S_Test)->AOV_Richness
summary(AOV_Richness)
write.csv(as.data.frame.summary.aovlist(summary(AOV_Richness)),"results/richness_AOV.csv")

library(rstatix)
test_s <- S_Test %>%
  group_by(Biome) %>%
  t_test(Richness ~ Origin) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_s 
write.csv(test_s,"results/richness_Tukey_biome.csv")
test_s <- test_s %>% add_xy_position(x = "Origin",step.increase = 1/14,fun="mean_ci",scales = "free")

test_s2 <- S_Test %>%
  group_by(Origin) %>%
  t_test(Richness ~ Biome) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_s2 
write.csv(test_s2,"results/richness_Tukey_origin.csv")

test_s[test_s$group1=="Exotic",]$xmin<-2
test_s[test_s$group2=="Exotic",]$xmax<-2
test_s[test_s$group1=="Native",]$xmin<-1
test_s[test_s$group2=="Native",]$xmax<-1

rm(AOV_Richness)
gc()

aggregate(S_Test$Richness,list(S_Test$Biome,S_Test$Mun,S_Test$State),mean)->All3
rbind(data.frame(Biome=All3$Group.1,Mun=All3$Group.2,State=All3$Group.3,Origin="All",Richness=All3$x),S_Test)->S_Test

aggregate(S_Test$Richness,list(S_Test$Biome,S_Test$Origin),FUN = function(x) c(mean = mean(x,na.rm=T), se = std.error(x,na.rm=T)))->MeanzRichness
MeanzRichness$lower<-MeanzRichness$x[,1]-(MeanzRichness$x[,2]*2)
MeanzRichness$upper<-MeanzRichness$x[,1]+(MeanzRichness$x[,2]*2)
colnames(MeanzRichness)[1:2]=c("Biome","Origin")

tiff("raw_figures/RichnessBars.tiff", units="in", width=7, height=3, res=900) #Save plot
{RichnessPlot=ggplot(MeanzRichness, aes(x=factor(Origin,levels=(c("Native","Exotic","All"))), y= x[,1], color= factor(Origin,levels=(c("Native","Exotic","All")))))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0,linewidth=.6,position=position_dodge(width=0.75)) +
  geom_point(size=1,position=position_dodge(width=0.75))+
  ggh4x::facet_grid2(.~Biome,scales="free",independent=T)+
  scale_x_discrete(labels=c("Native","Exotic","All"),name="Crop origin")+
  scale_color_manual(values=(c("#36AE7C","#EB5353","grey55")),name="Crop origin",labels=(c("Natives","Exotics","All crops")))+
  scale_y_continuous(name= "Crop richness",expand = expansion(mult = c(0.1, 0.1)))+
  stat_pvalue_manual(test_s, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)+
  theme_bw()+
    theme(legend.spacing.y = unit(0, 'cm'),strip.text = element_text(size=10,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(angle=90,family = "sans", colour = "black", size=12,vjust=0.5,hjust=1), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); RichnessPlot
}
dev.off()

##### Differences in PolDep ----
melt(BiomesData)->BiomesTest
BiomesTest[!grepl("All",BiomesTest$variable),]->BiomesTest
BiomesTest[grepl("PolDep",BiomesTest$variable),]->PolDep_Test
gsub("PolDep_Nat","Native",PolDep_Test$variable)->PolDep_Test$variable
gsub("PolDep_Exo","Exotic",PolDep_Test$variable)->PolDep_Test$variable
colnames(PolDep_Test)[4:5]=c("Origin","PolDep")
gsub("Amazônia","Amazon",PolDep_Test$Biome)->PolDep_Test$Biome
gsub("Pampa","Pampas",PolDep_Test$Biome)->PolDep_Test$Biome
gsub("Mata Atlântica","Atl. forest",PolDep_Test$Biome)->PolDep_Test$Biome

aov(PolDep~Biome*Origin,data = PolDep_Test)->AOV_PolDep
summary(AOV_PolDep)
write.csv(as.data.frame.summary.aovlist(summary(AOV_PolDep)),"results/poldep_AOV.csv")

library(rstatix)
test_poldep <- PolDep_Test %>%
  group_by(Biome) %>%
  t_test(PolDep ~ Origin) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_poldep 
write.csv(test_poldep,"results/poldep_Tukey_biome.csv")
test_poldep <- test_poldep %>% add_xy_position(x = "Origin",step.increase = 1/14,fun="mean_ci",scales = "free")

test_poldep2 <- PolDep_Test %>%
  group_by(Origin) %>%
  t_test(PolDep ~ Biome) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_poldep2 
write.csv(test_poldep2,"results/poldep_Tukey_origin.csv")

test_poldep[test_poldep$group1=="Exotic",]$xmin<-2
test_poldep[test_poldep$group2=="Exotic",]$xmax<-2
test_poldep[test_poldep$group1=="Native",]$xmin<-1
test_poldep[test_poldep$group2=="Native",]$xmax<-1
test_poldep[test_poldep$group1=="Exotic_NoSoy",]$xmin<-3
test_poldep[test_poldep$group2=="Exotic_NoSoy",]$xmax<-3

rm(AOV_PolDep)
gc()

aggregate(PolDep_Test$PolDep,list(PolDep_Test$Biome,PolDep_Test$Mun,PolDep_Test$State),mean)->All
rbind(data.frame(Biome=All$Group.1,Mun=All$Group.2,State=All$Group.3,Origin="All",PolDep=All$x),PolDep_Test)->PolDep_Test

aggregate(PolDep_Test$PolDep,list(PolDep_Test$Biome,PolDep_Test$Origin),FUN = function(x) c(mean = mean(x,na.rm=T), se = std.error(x,na.rm=T)))->MeanzPolDep
MeanzPolDep$lower<-MeanzPolDep$x[,1]-(MeanzPolDep$x[,2]*2)
MeanzPolDep$upper<-MeanzPolDep$x[,1]+(MeanzPolDep$x[,2]*2)
colnames(MeanzPolDep)[1:2]=c("Biome","Origin")

tiff("raw_figures/PolDepBars.tiff", units="in", width=10, height=3, res=900) #Save plot
{PolDep=ggplot(MeanzPolDep, aes(x=factor(Origin,levels=(c("Native","Exotic","Exotic_NoSoy","All"))), y= x[,1], color= factor(Origin,levels=(c("Native","Exotic","Exotic_NoSoy","All")))))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0,linewidth=.6,position=position_dodge(width=0.75)) +
  geom_point(size=1,position=position_dodge(width=0.75))+
    ggh4x::facet_grid2(.~Biome,scales="free",independent=T)+
    stat_pvalue_manual(test_poldep, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)+
    scale_x_discrete(labels=(c("Native","Exotic","No soy.","All")),name="Crop origin")+
  scale_color_manual(values=(c("#36AE7C","#bc4242","#EB5353","grey55")),name="Crop origin",labels=(c("Natives","Exotics","Exotics without\nsoybean","All crops")))+
  scale_y_continuous(name= "Proportion of production\ndependent on pollinators",expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(legend.spacing.y = unit(0, 'cm'),strip.text = element_text(size=10,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(angle=90,family = "sans", colour = "black", size=12,vjust=0.5,hjust=1), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); PolDep
}
dev.off()

#tiff("raw_figures/PolDepMaps.tiff", units="in", width=10, height=5, res=900) #Save plot
#ggarrange(PolDepMap_Nat,PolDepMap_Exo,PolDepMap_Exo_NoSoy,common.legend=TRUE,legend = "right",nrow = 1,ncol = 2)
#dev.off()

## Differences in PolRisk ----

BiomesTest[grepl("PolRisk",BiomesTest$variable),]->PolRisk_Test
gsub("PolRisk_Nat","Native",PolRisk_Test$variable)->PolRisk_Test$variable
gsub("PolRisk_Exo","Exotic",PolRisk_Test$variable)->PolRisk_Test$variable
colnames(PolRisk_Test)[4:5]=c("Origin","PolRisk")
gsub("Amazônia","Amazon",PolRisk_Test$Biome)->PolRisk_Test$Biome
gsub("Pampa","Pampas",PolRisk_Test$Biome)->PolRisk_Test$Biome
gsub("Mata Atlântica","Atl. forest",PolRisk_Test$Biome)->PolRisk_Test$Biome

aov(PolRisk~Biome*Origin,data = PolRisk_Test)->AOV_PolRisk
summary(AOV_PolRisk)
write.csv(as.data.frame.summary.aovlist(summary(AOV_PolRisk)),"results/pol_risk_AOV.csv")

library(rstatix)
test_polrisk <- PolRisk_Test %>%
  group_by(Biome) %>%
  t_test(PolRisk ~ Origin) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_polrisk 
write.csv(test_polrisk,"results/polrisk_Tukey_biome.csv")
test_polrisk <- test_polrisk %>% add_xy_position(x = "Origin",step.increase = 1/14,fun="mean_ci",scales = "free")

test_polrisk2 <- PolRisk_Test %>%
  group_by(Origin) %>%
  t_test(PolRisk ~ Biome) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_polrisk2 
write.csv(test_polrisk2,"results/polrisk_Tukey_origin.csv")

test_polrisk[test_polrisk$group1=="Exotic",]$xmin<-2
test_polrisk[test_polrisk$group2=="Exotic",]$xmax<-2
test_polrisk[test_polrisk$group1=="Native",]$xmin<-1
test_polrisk[test_polrisk$group2=="Native",]$xmax<-1
test_polrisk[test_polrisk$group1=="Exotic_NoSoy",]$xmin<-3
test_polrisk[test_polrisk$group2=="Exotic_NoSoy",]$xmax<-3

aggregate(PolRisk_Test$PolRisk,list(PolRisk_Test$Biome,PolRisk_Test$Mun,PolRisk_Test$State),mean)->All2
rbind(data.frame(Biome=All2$Group.1,Mun=All2$Group.2,State=All2$Group.3,Origin="All",PolRisk=All2$x),PolRisk_Test)->PolRisk_Test

aggregate(PolRisk_Test$PolRisk,list(PolRisk_Test$Biome,PolRisk_Test$Origin),FUN = function(x) c(mean = mean(x,na.rm=T), se = std.error(x,na.rm=T)))->MeanzPolRisk
MeanzPolRisk$lower<-MeanzPolRisk$x[,1]-(MeanzPolRisk$x[,2]*2)
MeanzPolRisk$upper<-MeanzPolRisk$x[,1]+(MeanzPolRisk$x[,2]*2)
colnames(MeanzPolRisk)[1:2]=c("Biome","Origin")

tiff("raw_figures/PolRiskBars.tiff", units="in", width=9.8, height=3, res=900) #Save plot
{PolRisk=ggplot(MeanzPolRisk, aes(x=factor(Origin,levels=(c("Native","Exotic","Exotic_NoSoy","All"))), y= x[,1], color= factor(Origin,levels=(c("Native","Exotic","Exotic_NoSoy","All")))))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0,linewidth=.6,position=position_dodge(width=0.75)) +
  geom_point(size=1,position=position_dodge(width=0.75))+
  ggh4x::facet_grid2(.~Biome,scales="free",independent=T)+
  scale_x_discrete(labels=c("Native","Exotic","No soy.","All"),name="Crop origin")+
  stat_pvalue_manual(test_polrisk, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)+
  scale_color_manual(values=(c("#36AE7C","#bc4242","#EB5353","grey55")),name="Crop origin",labels=(c("Natives","Exotics","Exotics without\nsoybean","All crops")))+
  scale_y_continuous(name= "Pollinator shortage risk",expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(strip.text = element_text(size=10,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); PolRisk
}
dev.off()

#tiff("raw_figures/PolRiskMaps.tiff", units="in", width=10, height=5, res=900) #Save plot
#ggarrange(PolRiskMap_Nat,PolRiskMap_Exo,common.legend=TRUE,legend = "right",nrow = 1,ncol = 2)
#dev.off()

#tiff("raw_figures/Figure 4.tiff", units="in", width=9, height=3, res=900) #Save plot
#ggarrange(PolDep,PolRisk,common.legend = TRUE,legend = "right",label.x = 0.27,label.y = 0.98,nrow = 1,ncol = 2,hjust = 1)
#dev.off()

###### Crop area -----
na.omit(AllMuns_NoZero)->AllMuns_NoZero

sum(AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Nat",]$AreaPlt)*100/sum(AllMuns_NoZero$AreaPlt) #18.43% Native
sum(AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Exo",]$AreaPlt)*100/sum(AllMuns_NoZero$AreaPlt) #81.57% Exotic
sum(AllMuns_NoZero[grepl("Soja",AllMuns_NoZero$Crop),]$AreaPlt)*100/sum(AllMuns_NoZero$AreaPlt) #46.04% Soybean
NoSoybean<-AllMuns_NoZero[!grepl("Soja",AllMuns_NoZero$Crop),]
sum(NoSoybean[NoSoybean$crop_origin=="Nat",]$AreaPlt)*100/sum(NoSoybean$AreaPlt)
sum(NoSoybean[NoSoybean$crop_origin=="Exo",]$AreaPlt)*100/sum(NoSoybean$AreaPlt)

length(unique(AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Nat",]$Crop))
length(unique(AllMuns_NoZero[AllMuns_NoZero$crop_origin=="Exo",]$Crop))

nrow(pol_depend[!is.na(pol_depend$pol_dep),])/nrow(pol_depend)
nrow(pol_depend[is.na(pol_depend$pol_dep),])
nrow(pol_depend[is.na(pol_depend$pol_dep),][pol_depend[is.na(pol_depend$pol_dep),]$origin=="Nat",])/nrow(pol_depend[is.na(pol_depend$pol_dep),])
pol_depend[!is.na(pol_depend$pol_dep),]->only_avai
nrow(only_avai[!only_avai$pol_dep=="no increase",])/nrow(only_avai)
pol_depend[pol_depend$origin=="Nat",][!is.na(pol_depend[pol_depend$origin=="Nat",]$pol_dep),]->NatPol
nrow(NatPol[!NatPol$pol_dep=="no increase",])/nrow(NatPol)
nrow(NatPol[NatPol$pol_dep=="essential",])/nrow(NatPol)
nrow(NatPol[NatPol$pol_dep=="great",])/nrow(NatPol)
nrow(NatPol[NatPol$pol_dep=="modest",])/nrow(NatPol)
nrow(NatPol[NatPol$pol_dep=="little",])/nrow(NatPol)
pol_depend[pol_depend$origin=="Exo",][!is.na(pol_depend[pol_depend$origin=="Exo",]$pol_dep),]->ExoPol
nrow(ExoPol[!ExoPol$pol_dep=="no increase",])/nrow(ExoPol)
nrow(ExoPol[ExoPol$pol_dep=="essential",])/nrow(ExoPol)
nrow(ExoPol[ExoPol$pol_dep=="great",])/nrow(ExoPol)
nrow(ExoPol[ExoPol$pol_dep=="modest",])/nrow(ExoPol)
nrow(ExoPol[ExoPol$pol_dep=="little",])/nrow(ExoPol)

melt(table(pol_depend$origin,pol_depend$pol_dep))->PolPerc

PolPerc_Plot= ggplot(PolPerc, aes(y=value, x=Var1, fill= factor(Var2,levels=c("no increase","little","modest","great","essential"))))+
  geom_bar(position= "fill", stat= "identity")+
  coord_flip()+
  scale_fill_manual(labels=c("None","Little","Modest","Great","Essential"),values = moma.colors("VanGogh", 5), name= "Pollinator\ndependence")+
  scale_x_discrete(name= "Plant origin", labels= c("Exotic","Native"))+
  scale_y_continuous(name= "Proportion of crop species (%)")+
  theme_bw()+
  theme(legend.position = "right",legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=14), legend.text = element_text(family ="sans", size=12), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14));PolPerc_Plot

tiff("raw_figures/Figure 3.tiff", units="in", width=6, height=4, res=900)
ggarrange(sistreplot,PolPerc_Plot,ncol = 1,nrow = 2,labels = c("(a)","(b)"),hjust = 0,label.x = 0.09,label.y = 0.97,align = "hv")
dev.off()

BiomesTest[grepl("Area_",BiomesTest$variable),]->Area_Test
gsub("Area_Nat","Native",Area_Test$variable)->Area_Test$variable
gsub("^Area_Exo$","Exotic",Area_Test$variable)->Area_Test$variable
colnames(Area_Test)[4:5]=c("Origin","Area")
gsub("Amazônia","Amazon",Area_Test$Biome)->Area_Test$Biome
gsub("Pampa","Pampas",Area_Test$Biome)->Area_Test$Biome
gsub("Mata Atlântica","Atl. forest",Area_Test$Biome)->Area_Test$Biome

aov(Area~Biome*Origin,data = Area_Test)->AOV_Area
summary(AOV_Area)
write.csv(as.data.frame.summary.aovlist(summary(AOV_Area)),"results/crop_area_AOV.csv")

library(rstatix)
test_area <- Area_Test %>%
  group_by(Biome) %>%
  t_test(Area ~ Origin) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_area 
write.csv(test_area,"results/crop_area_Tukey_biome.csv")
test_area <- test_area %>% add_xy_position(x = "Origin",step.increase = 1/14,fun="mean_ci",scales = "free")

test_area2 <- Area_Test %>%
  group_by(Origin) %>%
  t_test(Area ~ Biome) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
test_area2 
write.csv(test_area2,"results/crop_area_Tukey_origin.csv")

test_area[test_area$group1=="Exotic",]$xmin<-2
test_area[test_area$group2=="Exotic",]$xmax<-2
test_area[test_area$group1=="Native",]$xmin<-1
test_area[test_area$group2=="Native",]$xmax<-1
test_area[test_area$group1=="Area_ExoNoSoy",]$xmin<-3
test_area[test_area$group2=="Area_ExoNoSoy",]$xmax<-3

aggregate(Area_Test$Area,list(Area_Test$Biome,Area_Test$Mun,Area_Test$State),mean)->All4
rbind(data.frame(Biome=All4$Group.1,Mun=All4$Group.2,State=All4$Group.3,Origin="All",Area=All4$x),Area_Test)->Area_Test

aggregate(Area_Test$Area,list(Area_Test$Biome,Area_Test$Origin),FUN = function(x) c(mean = mean(x,na.rm=T), se = std.error(x,na.rm=T)))->MeanzArea
MeanzArea$lower<-MeanzArea$x[,1]-(MeanzArea$x[,2]*2)
MeanzArea$upper<-MeanzArea$x[,1]+(MeanzArea$x[,2]*2)
colnames(MeanzArea)[1:2]=c("Biome","Origin")

tiff("raw_figures/Area_Bars.tiff", units="in", width=9, height=3, res=900) #Save plot
ErrorPlots=ggplot(MeanzArea, aes(x=factor(Origin,levels=c("Native","Exotic","Area_ExoNoSoy","All")), y= x[,1], color= factor(Origin,levels=c("Native","Exotic","Area_ExoNoSoy","All"))))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0,linewidth=.6,position=position_dodge(width=0.75)) +
  geom_point(size=2,position=position_dodge(width=0.75))+
  ggh4x::facet_grid2(.~Biome,scales="free",independent=T)+
  stat_pvalue_manual(test_area, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)+
  scale_x_discrete(labels=c("Native","Exotic","No soy.","All"),name="Crop origin")+
  scale_color_manual(values=c("#36AE7C","#bc4242","#EB5353","grey55"),name="Crop origin",labels=c("Native","Exotic","Exotic without\nsoybean","All crops"))+
  scale_y_continuous(name= "Crop area (hectare)",expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(strip.text = element_text(size=10,family="sans",face = "bold"),strip.background = element_rect(fill="white"), legend.text.align = 0,legend.title.align=0,legend.position = "none", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA")); ErrorPlots
dev.off()

aggregate(AllMuns_Nat$AreaPlt,list(AllMuns_Nat$Município),sum,na.rm=T)->AreaPlt_Nat
aggregate(AllMuns_Exo$AreaPlt,list(AllMuns_Exo$Município),sum,na.rm=T)->AreaPlt_Exo
aggregate(AllMuns_Exo[!grepl("Soja",AllMuns_Exo$Crop),]$AreaPlt,list(AllMuns_Exo[!grepl("Soja",AllMuns_Exo$Crop),]$Município),sum,na.rm=T)->AreaPlt_Exo_NoSoy
aggregate(AllMuns$AreaPlt,list(AllMuns$Município),sum,na.rm=T)->AreaPlt_All
colnames(AreaPlt_All)=c("Mun","AreaPlt")
AreaPlt_All$AreaNative<-NA
AreaPlt_All$AreaExotic<-NA
AreaPlt_All$AreaExotic_NoSoy<-NA

AreaPlt_All[match(AreaPlt_Nat$Group.1,AreaPlt_All$Mun),]$AreaNative<-AreaPlt_Nat$x
AreaPlt_All[match(AreaPlt_Exo$Group.1,AreaPlt_All$Mun),]$AreaExotic<-AreaPlt_Exo$x
AreaPlt_All[match(AreaPlt_Exo_NoSoy$Group.1,AreaPlt_All$Mun),]$AreaExotic_NoSoy<-AreaPlt_Exo_NoSoy$x
AreaPlt_All$AreaNative[is.na(AreaPlt_All$AreaNative)]<-0
AreaPlt_All$AreaExotic[is.na(AreaPlt_All$AreaExotic)]<-0

AreaPlt_All$PercNat<-AreaPlt_All$AreaNative/AreaPlt_All$AreaPlt
AreaPlt_All$PercExo<-AreaPlt_All$AreaExotic/AreaPlt_All$AreaPlt
AreaPlt_All$PercExo_NoSoy<-AreaPlt_All$AreaExotic_NoSoy/AreaPlt_All$AreaPlt

BRMun$PercAreaPlt_Nat<-NA
BRMun$PercAreaPlt_Exo<-NA
BRMun$PercAreaPlt_Exo_NoSoy<-NA
BRMun$AreaPlt<-NA
BRMun[match(AreaPlt_All$Mun,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$PercAreaPlt_Nat<-AreaPlt_All$PercNat
BRMun[match(AreaPlt_All$Mun,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$PercAreaPlt_Exo<-AreaPlt_All$PercExo
BRMun[match(AreaPlt_All$Mun,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$PercAreaPlt_Exo_NoSoy<-AreaPlt_All$PercExo_NoSoy
BRMun[match(AreaPlt_All$Mun,paste0(BRMun$NM_MUN," (",BRMun$SIGLA_UF,")")),]$AreaPlt<-AreaPlt_All$AreaPlt

AreaPltMap1<-ggplot() +
  geom_sf(data = BRMun,aes(color=PercAreaPlt_Exo,fill=PercAreaPlt_Exo)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_color_gradient(low = "green",high = "red",name="Crop prevalence",breaks=c(0,0.5,1),labels=c("More native",0.5,"More exotic"))+
  scale_fill_gradient(low = "green",high = "red",name="Crop prevalence",breaks=c(0,0.5,1),labels=c("More native",0.5,"More exotic"))+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))

AreaPltMap2<-ggplot() +
  geom_sf(data = BRMun,aes(color=PercAreaPlt_Exo_NoSoy,fill=PercAreaPlt_Exo_NoSoy)) +
  geom_sf(data=Biomes,fill=NA,color="black",linewidth=0.2)+
  scale_x_continuous(limits=c(-73.99045,-35.0000))+
  scale_color_gradient(low = "green",high = "red",name="Crop prevalence",breaks=c(0,0.5,1),labels=c("More native",0.5,"More exotic"))+
  scale_fill_gradient(low = "green",high = "red",name="Crop prevalence",breaks=c(0,0.5,1),labels=c("More native",0.5,"More exotic"))+
  theme_bw()+
  theme(legend.text.align = 0,legend.position = "right", legend.background = element_rect(fill="NA", colour = "NA"), legend.title = element_text(face="bold", family = "sans", colour = "black", size=12), legend.text = element_text(family ="sans", size=10), axis.text.x = element_text(family = "sans", colour = "black", size=12), axis.text.y = element_text(family = "sans", colour = "black", size=12), axis.title = element_text(face="bold",family = "sans", size = 14), panel.border = element_rect(colour = "black", fill = "NA"))


tiff("raw_figures/CropAreaMaps.tiff", units="in", width=10, height=5, res=900) #Save plot
ggarrange(AreaPltMap1,AreaPltMap2,nrow = 1,ncol = 2,common.legend = T,legend = "right")
dev.off()
