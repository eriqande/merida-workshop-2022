---
output:
  word_document: default
  html_document: default
---
#R Code for Echeverria Caro A, Feldman RE, Bahn V.
#Geographic context is a stronger driver of bird diversity during migration than primary productivity.

##Required packages
```{r eval=TRUE}
library(tidyverse)
library(viridis)
library(ape)
library(modelr)
library(ggpubr)
```

##Import bird data.
###Data available on the Dryad Digital Repository FORTHCOMING

data.all<-read_csv("data.all.csv"))

###Mean center and standardize all predictor variables and calculate max species richness per eBird checklist

```{r eval=FALSE}
data.std<-data.all.2%>%select(TEMPORADA, LATITUDE, LONGITUDE, DURATION.MINUTES, OBSERVATION.DATE, PERIODO, SCIENTIFIC.NAME, 
                            EFFORT.DISTANCE.KM, Dist_km_new, NDVI.med, EVI.med, HET.med)%>%
  group_by(TEMPORADA, LATITUDE, LONGITUDE, DURATION.MINUTES, OBSERVATION.DATE, PERIODO, 
           EFFORT.DISTANCE.KM, Dist_km_new, NDVI.med, EVI.med, HET.med)%>%
  summarize(Riqueza=n_distinct(SCIENTIFIC.NAME))%>%
  pivot_longer(cols=c(NDVI.med, EVI.med, HET.med, Dist_km_new, EFFORT.DISTANCE.KM), names_to="Index", values_to="Value")%>%
  filter(Value!=0)%>%group_by(TEMPORADA, Index)%>%nest()%>%mutate(Std=map(data, function(x){as.numeric(scale(x$Value))}))%>%
  unnest(cols=c(data, Std))%>%select(-Value)%>%
  pivot_wider(names_from=Index, values_from=Std)%>%group_by(TEMPORADA, LATITUDE, LONGITUDE, PERIODO)%>%summarize(Riqueza.m=max(Riqueza),
                                                                                                                 NDVI.m=mean(NDVI.med),
                                                                                                                 EVI.m=mean(EVI.med),
                                                                                                                 HET.m=mean(HET.med),
                                                                                                                 DIST.m=mean(Dist_km_new),
                                                                                                                 EFFORT.m=mean(EFFORT.DISTANCE.KM))%>%
  pivot_longer(cols=c(NDVI.m, EVI.m, HET.m), names_to="Index", values_to="Value")%>%mutate(Riqueza.m=as.integer(floor(Riqueza.m)))

```

###Function to calculate pseudo-R2 https://stackoverflow.com/questions/57319130/purrrmap-and-glm-issues-with-call

```{r eval=FALSE}
get.r2<-function(x){
  
  L.base<-
    logLik(
      glm(formula = reformulate('1',"Riqueza.m"),
          data=x$data,
          family = x$family))
  
  n<-length(x$residuals)
  
  L.full<-logLik(x)
  D.full <- -2 * L.full
  D.base <- -2 * L.base
  return(data.frame(Nagelkerke = (1 - exp((D.full - D.base)/n))/(1 - exp(-D.base/n))))
  
}
```

###Function to run 5-fold cross validation and calculate root mean square error 

```{r eval=FALSE}
get.cv<-function(x){
  x%>%crossv_kfold(k = 5)%>%
    mutate(model = map(train, ~glm(unique(x$Form), family=poisson, data=.)),
           rmse_all_models = map2_dbl(model, test, ~rmse(.x, .y)))%>%pull(rmse_all_models)%>%mean()
}
```

###Function to generate data frame used for predicting species richness from models. Predictions are made for a continuous range of primary productivity and distance to coast values

```{r eval=FALSE}
predict.frame<-function(x){
  expand_grid(Value=seq(floor(min(x$Value)), max(ceiling(x$Value)), 0.1), 
              DIST.m=seq(floor(min(x$DIST.m)), ceiling(max(x$DIST.m)), 0.1), 
              EFFORT.m=0, PERIODO=unique(x$PERIODO))
}

```

###Function to calculate spatial autocorrelation (Moran's I) for the residuals of each model

```{r eval=FALSE}
get.spatauto<-function(x,y){
w<-as.matrix(dist(cbind(x$LATITUDE, x$LATITUDE)))
w1<-1/w
diag(w1)<-0
w1[is.infinite(w1)]<-0
m<-Moran.I(y$residuals,w1)
tibble(Observed=m[[1]], Expected=m[[2]], SD=m[[3]], P=m[[4]])
}
```

###Set-up the data frame to run full and reduced models

```{r eval=FALSE}
data.models<-data.std%>%mutate(Mod1="Riqueza.m~EFFORT.m+DIST.m+Value+PERIODO+DIST.m:PERIODO+Value:PERIODO",
                                Mod2="Riqueza.m~EFFORT.m+DIST.m+Value+PERIODO+DIST.m:PERIODO",
                                Mod3="Riqueza.m~EFFORT.m+DIST.m+Value+PERIODO+Value:PERIODO",
                                Mod4="Riqueza.m~EFFORT.m+DIST.m+Value+PERIODO",
                                Mod5="Riqueza.m~EFFORT.m+DIST.m+PERIODO+DIST.m:PERIODO",
                                Mod6="Riqueza.m~EFFORT.m+Value+PERIODO+Value:PERIODO",
                                Mod7="Riqueza.m~EFFORT.m+DIST.m+PERIODO",
                                Mod8="Riqueza.m~EFFORT.m+Value+PERIODO")%>%pivot_longer(cols=c(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8), names_to="Type", values_to="Form")%>%
  filter(!is.na(Value))
````

###Run models

```{r eval=FALSE}
data.results<-data.models%>%group_by(TEMPORADA, Index, Type)%>%nest()%>%
  mutate(Model=map(data, ~glm(unique(.$Form), family="poisson", data=.)))%>%
  mutate(Model.results=Model %>% map(tidy), Stats=map(Model, AIC))%>%mutate(Rsq=map(Model,get.r2), RMSE=map(data, get.cv))#%>%
  mutate(Moran.test=map2(data, Model, get.spatauto))
```

###Predict species richness for models 2 and models 5 (best models for fall and spring, respectively)

```{r eval=FALSE}
data.predictions<-data.results%>%mutate(Final.model=paste(TEMPORADA, Index, Type, sep="."))%>%
  filter(Final.model%in%c("OTO?O.NDVI.m.Mod2", "PRIMAVERA.NDVI.m.Mod5"))%>%select(Model, data)%>%
  mutate(Predict.data=map(data, predict.frame))%>%
  mutate(Predictions = map2(.x=Model, .y=Predict.data, .f=~augment(.x, newdata=.y, type.predict="response", se_fit=TRUE, interval="prediction")))

```

##Graph relationship with Distance and NDVI
###To display actual distance and NDVI values (unstandardized), need to calculate the mean and sd for the variables

```{r eval=FALSE}
data.means<-data.all.2%>%select(TEMPORADA, LATITUDE, LONGITUDE, DURATION.MINUTES, OBSERVATION.DATE, PERIODO, SCIENTIFIC.NAME, 
                              EFFORT.DISTANCE.KM, Dist_km_new, NDVI.med, EVI.med, HET.med)%>%
  group_by(TEMPORADA, LATITUDE, LONGITUDE, DURATION.MINUTES, OBSERVATION.DATE, PERIODO, 
           EFFORT.DISTANCE.KM, Dist_km_new, NDVI.med, EVI.med, HET.med)%>%
  summarize(Riqueza=n_distinct(SCIENTIFIC.NAME))%>%
  pivot_longer(cols=c(NDVI.med, EVI.med, HET.med, Dist_km_new, EFFORT.DISTANCE.KM), names_to="Index", values_to="Value")%>%
  filter(Value!=0)%>%group_by(TEMPORADA, Index)%>%
  summarize(Mean=mean(Value), SD=sd(Value))%>%filter(Index%in%c("Dist_km_new", "NDVI.med"))%>%
  pivot_wider(names_from="Index", values_from=c(Mean, SD))
```

###Join with original data to display as points on the graphs

```{r eval=FALSE}
data.raw<-data.all.2%>%select(TEMPORADA, LATITUDE, LONGITUDE, DURATION.MINUTES, OBSERVATION.DATE, PERIODO, SCIENTIFIC.NAME, 
                                EFFORT.DISTANCE.KM, Dist_km_new, NDVI.med, EVI.med, HET.med)%>%
  group_by(TEMPORADA, LATITUDE, LONGITUDE, DURATION.MINUTES, OBSERVATION.DATE, PERIODO, 
           EFFORT.DISTANCE.KM, Dist_km_new, NDVI.med, EVI.med, HET.med)%>%
  summarize(Riqueza=n_distinct(SCIENTIFIC.NAME))%>%rename(.fitted=Riqueza, NDVI=NDVI.med, Distance=Dist_km_new)%>%
  mutate(Index=NA, Type=NA, Value=NA, DIST.m=NA, .se.fit=NA)%>%ungroup()%>%select(TEMPORADA, Index, Type, Value, DIST.m, PERIODO, .fitted, .se.fit, NDVI, Distance)

```

###Prepare data for graphing

```{r eval=FALSE}
graphdata<-data.predictions%>%select(Predictions)%>%unnest(Predictions)%>%full_join(., data.means, by="TEMPORADA")%>%
  mutate(NDVI=Value*SD_NDVI.med+Mean_NDVI.med, 
         Distance=DIST.m*SD_Dist_km_new+Mean_Dist_km_new)%>%select(TEMPORADA, Index, Type, Value, DIST.m, PERIODO, .fitted, .se.fit, NDVI, Distance)%>%
  bind_rows(., data.raw)
```

##Figure 1. Relationship between NDVI and species richness for fall migration

```{r eval=FALSE}
ndvi.fall<-ggplot(graphdata%>%filter(TEMPORADA=="OTO?O", DIST.m==0), 
                      aes(x=NDVI, y=.fitted, col=factor(PERIODO))) +  geom_line(size=2)
ndvi.fall<-ndvi.fall+geom_point(data=graphdata%>%filter(TEMPORADA=="OTO?O", is.na(Value), .fitted!=26),#Remove outlier to facilitate visualizing 
                                        aes(x=NDVI, y=.fitted, col=factor(PERIODO)), size=1)
ndvi.fall<-ndvi.fall+scale_y_continuous(name="Species richness\n")
ndvi.fall<-ndvi.fall+scale_x_continuous(name="\nNDVI")
ndvi.fall<-ndvi.fall+scale_color_viridis_d(name="",labels = c("29/Aug-13/Sep", "14/Sep-29/Sep","30/Sep-15/Oct","16/Oct-31/Oct","01/Nov-16/Nov","17/Nov-02/Dec"),aesthetics=c("color", "fill"))
ndvi.fall<-ndvi.fall+geom_ribbon(aes(ymin=.fitted-.se.fit, ymax=.fitted+.se.fit, fill=factor(PERIODO)), 
                                         show.legend=FALSE, color=NA, alpha=0.3)
ndvi.fall<-ndvi.fall+theme_classic()+ theme(axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"), 
                                            text = element_text(size=12))

ndvi.fall
ggsave("Figure.1.pdf", device="pdf", width=6, height=3, units="in", dpi=600)
```

##Figure 2. Relationship between distance to the coast and species richness for fall and spring migration

```{r eval=FALSE}
distance.fall<-ggplot(graphdata%>%filter(TEMPORADA=="OTO?O", Value==0), 
                      aes(x=Distance, y=.fitted, col=factor(PERIODO))) +  geom_line(size=2)
distance.fall<-distance.fall+geom_point(data=graphdata%>%filter(TEMPORADA=="OTO?O", is.na(Value), .fitted!=26),#Remove outlier to facilitate visualizing 
                      aes(x=Distance, y=.fitted, col=factor(PERIODO)), size=1)
distance.fall<-distance.fall+scale_y_continuous(name="Species richness\n")#, 
                                                #breaks = c(2,4,6,8,10))
distance.fall<-distance.fall+scale_x_continuous(name="\nDistance to coast (km)")

distance.fall<-distance.fall+scale_color_viridis_d(name="",labels = c("29/Aug-13/Sep", "14/Sep-29/Sep","30/Sep-15/Oct","16/Oct-31/Oct","01/Nov-16/Nov","17/Nov-02/Dec"),aesthetics=c("color", "fill"))

distance.fall<-distance.fall+geom_ribbon(aes(ymin=.fitted-.se.fit, ymax=.fitted+.se.fit, fill=factor(PERIODO)), 
                                         show.legend=FALSE, color=NA, alpha=0.3)
distance.fall<-distance.fall+theme_classic()+ theme(axis.line.x=element_blank(), axis.line.y=element_line(colour="black"), 
                                                     text = element_text(size=12), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

distance.fall
#ggsave("Distance.fall.png", device="png", width=23.4, height=17.4, units="cm", dpi=300)


distance.spring<-ggplot(graphdata%>%filter(TEMPORADA=="PRIMAVERA", Value==0), 
                      aes(x=Distance, y=.fitted, col=factor(PERIODO))) +  geom_line(size=2)
distance.spring<-distance.spring+geom_point(data=graphdata%>%filter(TEMPORADA=="PRIMAVERA",is.na(Value)),#Remove outlier to facilitate visualizing 
                                        aes(x=Distance, y=.fitted, col=factor(PERIODO)), size=1)
distance.spring<-distance.spring+scale_y_continuous(name="Species richness\n")

distance.spring<-distance.spring+scale_x_continuous(name="\nDistance to coast (km)")

distance.spring<-distance.spring+scale_color_viridis_d(name="",labels = c("18/Feb-05/Mar", "06/Mar-21/Mar","22/Mar-06/Apr","07/Apr-22/Apr","23/Apr-08/May","09/May-24/May", "25/May-09/Jun"),aesthetics=c("color", "fill"))

distance.spring<-distance.spring+geom_ribbon(aes(ymin=.fitted-.se.fit, ymax=.fitted+.se.fit, fill=factor(PERIODO)), 
                                         show.legend=FALSE, color=NA, alpha=0.3)
distance.spring<-distance.spring+theme_classic()+ theme(axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"), 
                                                    text = element_text(size=12), axis.title.y=element_blank())

distance.spring
#ggsave("Distance.spring.png", device="png", width=23.4, height=17.4, units="cm", dpi=300)

##One figure combining both seasons
distance.both<-ggarrange(distance.fall, distance.spring, ncol = 1, labels=c("a", "b"), label.x=0.05,  
                               font.label=list(size=10), legend="right", align="v", common.legend=FALSE)
annotate_figure(distance.both, left="Species richness\n")
ggsave("Figure.2.pdf", device="pdf", width=6, height=6, units="in", dpi=600)
```

