####load packages####
library(spatstat)
library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
library(tibble)
library(MuMIn)


####load data####
ljcdt<-read.csv("./Lan5.28_analysisdata.csv") %>%
  mutate(.,sp.x=ifelse(sp.x=="白?","白臼",sp.x))
ljcenv<-read.csv("./欖仁溪環境資料.csv") %>%
  mutate(.,x1=(X-26)*10,y1=(Y-17)*10)
# splist<-read.csv("./TL50splist.csv") %>%
#   mutate(.,sp=ifelse(sp=="白?","白臼",sp))

####15 splist####
splist<-ljcdt %>%
  group_by(.,sp.x) %>%
  summarise(.,dead19=sum(dead18),surviving19=sum(live13==1 & live18==1)) %>%
  ungroup(.) %>%
  filter(.,dead19>=15 & surviving19>=15)

####envDT####
elevdt<-ljcenv %>%
  filter(.,X<=49) %>%
  select(.,x1,y1,Eleva)
colnames(elevdt)[1]<-"x"
colnames(elevdt)[2]<-"y"
elev<-as.im(elevdt)
unitname(elev) <- c("metre", "metres")

angdt<-ljcenv %>%
  filter(.,X<=49) %>%
  select(.,x1,y1,Ang.90)
colnames(angdt)[1]<-"x"
colnames(angdt)[2]<-"y"
ang90<-as.im(angdt)
unitname(ang90) <- c("metre", "metres")

aspdt<-ljcenv %>%
  filter(.,X<=49) %>%
  select(.,x1,y1,Asp.225)
colnames(aspdt)[1]<-"x"
colnames(aspdt)[2]<-"y"
asp225<-as.im(aspdt)
unitname(asp225) <- c("metre", "metres")

slopedt<-ljcenv %>%
  filter(.,X<=49) %>%
  select(.,x1,y1,Slope)
colnames(slopedt)[1]<-"x"
colnames(slopedt)[2]<-"y"
slope<-as.im(slopedt)
unitname(slope) <- c("metre", "metres")

convexdt<-ljcenv %>%
  filter(.,X<=49) %>%
  select(.,x1,y1,Convex)
colnames(convexdt)[1]<-"x"
colnames(convexdt)[2]<-"y"
convex<-as.im(convexdt)
unitname(convex) <- c("metre", "metres")

ljc.extra<-solist(elev=elev,ang90=ang90,asp225=asp225,slope=slope,convex=convex)

####spPattern####
spname<-"小葉赤楠"
ljcpoint<-ljcdt %>%
  filter(.,x4<=230,y4<=210,sp.x==spname,daed18==1) %>%
  select(.,x4,y4)
ljcpattern<-ppp(ljcpoint$x4,ljcpoint$y4,c(0,230),c(0,210))

####all point pcf####
ljcpoint<-ljcdt %>% filter(.,x4<=230,y4<=210,dead18==1) %>% select(.,x4,y4)
ljcpattern<-ppp(ljcpoint$x4,ljcpoint$y4,c(0,230),c(0,210))
E<-envelope(ljcpattern, pcf, nsim = 99, divisor="d", correction="Ripley")
png(filename = paste0("pcf-dead19-all.png"),width = 800,height = 625)
plot(E,legend=F,main = "")
dev.off()

####ppm AIC/BIC####
model_1<-ppm(ljcpattern, ~elev+ang90+asp225+slope+convex)
model_2<-ppm(ljcpattern, ~ang90+asp225)
model_3<-ppm(ljcpattern, ~ang90+asp225+convex)
model_4<-ppm(ljcpattern, ~ang90+asp225+slope+convex)
model_5<-ppm(ljcpattern, ~asp225)
model_6<-ppm(ljcpattern, ~ang90)
model_7<-ppm(ljcpattern, ~elev)
model_8<-ppm(ljcpattern, ~elev+slope+convex)
model_9<-ppm(ljcpattern, ~elev+ang90+slope+convex)
model_10<-ppm(ljcpattern, ~ang90*asp225)

AIC_seq<-model.sel(model_1,model_2,model_3,model_4,model_5,model_6,model_7,model_8,model_9,model_10,rank=AIC)
AIC_seq

BIC_seq<-model.sel(model_1,model_2,model_3,model_4,model_5,model_6,model_7,model_8,model_9,model_10,rank=BIC)
BIC_seq

####ppm(env)####
ljcpoint<-ljcdt %>% filter(.,x4<=230,y4<=210,dead18==1) %>% select(.,x4,y4)
ljcpattern<-ppp(ljcpoint$x4,ljcpoint$y4,c(0,230),c(0,210))
fit<-ppm(ljcpattern, ~elev+ang90+asp225+slope+convex,covariates = ljc.extra)
E<-envelope(ljcpattern, pcfinhom, nsim = 99,funargs = list(lambda=fit),correction="Ripley", divisor="d")
png(filename = paste0("pcfinhom-dead19-allenv.png"),width = 800,height = 625)
plot(E,legend=F,main = "")
dev.off()

####rlabel+env ppm####
surviving19df<-ljcdt %>%
  filter(.,x4<=230,y4<=210,live13==1 & live18==1) %>%
  select(.,x4,y4) %>%
  mutate(.,mark="A")
dead19df<-ljcdt %>%
  filter(.,x4<=230,y4<=210,dead18==1) %>%
  select(.,x4,y4) %>%
  mutate(.,mark="D")
ljcpatterndf<-bind_rows(surviving19df,dead19df)
ljcpattern<-ppp(ljcpatterndf$x4,ljcpatterndf$y4,c(0,230),c(0,210),marks = as.factor(ljcpatterndf$mark))
fit<-ppm(ljcpattern, ~elev+ang90+asp225+slope+convex,covariates = ljc.extra)
Epcfcrossinhom<-envelope(ljcpattern, pcfcross.inhom, nsim = 99,funargs = list(lambdaX=fit),correction="Ripley", divisor="d", i="D", j="D",simulate = expression(rlabel(ljcpattern)))
png(filename = paste0("pcfcrossinhom-dead19-allenv.png"),width = 800,height = 625)
plot(Epcfcrossinhom,legend=F,main = "")
dev.off()


####***dense+env (intensity) ppm***####

set.seed(1234)
Rscale<-30
ljcpoint<-ljcdt %>% filter(.,x4<=230,y4<=210,dead18==1) %>% select(.,x4,y4)
ljcpattern<-ppp(ljcpoint$x4,ljcpoint$y4,c(0,230),c(0,210))
intensityR<-density(ljcpattern,sigma = Rscale, kernel = "epanechnikov",positive = T)
# intensityR$v<-ifelse(intensity10$v<0,0,intensity10$v)
ljc.extra<-solist(elev=elev,ang90=ang90,asp225=asp225,slope=slope,convex=convex,intensityR=intensityR)
fit<-ppm(ljcpattern, ~elev+ang90+asp225+slope+convex+intensityR,covariates = ljc.extra)
intensitymap<-intensity(fit)
EpcfI<-envelope(ljcpattern, pcf, nsim = 99,correction="Ripley", divisor="d",simulate = expression(rpoispp(intensitymap)))
png(filename = paste0("pcf-dead19(intensity",Rscale,")-allenv.png"),width = 800,height = 625)
plot(EpcfI,legend=F,main = "")
dev.off()


####point pattern LOOP####
for (i in 1:nrow(splist)) {
  ljcpoint<-ljcdt %>% 
    filter(.,X1<49,Y1<38,sp.x==splist[i,1],dead18==1) %>%
    select(.,x4,y4)
  ljcpattern<-ppp(ljcpoint$x4,ljcpoint$y4,c(0,230),c(0,210))
  
  dense<-density(ljcpattern,sigma = 10)
  dense$v<-ifelse(dense$v<0,0,dense$v)
  fit<-ppm(ljcpattern, ~elev+ang90+asp225+slope+convex,covariates = ljc.extra)
  E<-envelope(ljcpattern, Linhom, nsim = 99,funargs = list(lambda=fit), simulate = expression(rpoispp(dense)),correction="best")
  png(filename = paste0("dead19-",i,"-",splist[i,1],".png"),width = 800,height = 625)
  plot(E,.-r~r,legend=F,main = splist[i,1])
  dev.off()
}
