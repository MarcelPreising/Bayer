if(!require(MASS)) install.packages("MASS")
require(MASS)
library(MASS)

if(!require(pscl)) install.packages("pscl")
require(pscl)
library(pscl)

if(!require(ggplot2)) install.packages("ggplot2")
require(ggplot2)
library(ggplot2)
require(gridExtra)

devtools::install_github("hrbrmstr/waffle")
library(waffle)

if(!require(rpart)) install.packages("rpart")
require(rpart)
library(rpart)


### Read in Data:
first <- read.csv("C:/.../data/randomization.csv",header=TRUE)
second <- read.csv("C:/.../data/subject.csv",header=TRUE)
third <- read.csv("C:/.../data/efficacy.csv",header=TRUE)

### Data in one Dataframe:
Data <- matrix(NA,nrow = dim(first[1]), ncol = dim(first)[2]+dim(second)[2]+dim(third)[2] -2)
Data <- as.data.frame(Data)
n <- nrow(Data)
K <- ncol(Data)-1

# Use order in table from scenario description:
Data[,1] <- first[,1]
Data[,2] <- first[,2]
Data[,3] <- second[,2]
Data[,4] <- second[,3]
Data[,5] <- second[,6]
Data[,6] <- second[,4]
Data[,7] <- second[,5]
Data[,8] <- third[,2]
Data[,9] <- third[,3]


### Names of the variables:
varnames <- c(names(first),names(second)[c(2,3,6,4,5)],names(third)[2:3])
# Column names:
colnames(Data) <- varnames

Data_Summary <- summary(Data)



### Wandel Treatment in Dummy-Variable um, 1=Active, 0=Placebo:
Data[,2] <- as.numeric(Data[,2])
for(i in 1:n){
  if(Data[i,2] == 2) Data[i,2] <- 0
}

### Wandel Tissue-Use in Dummy-Variable um, 1=High, 0=Medium:
Data[,6] <- as.numeric(Data[,6])
for(i in 1:n){
  if(Data[i,6] == 2) Data[i,6] <- 0
}


### Welche Individuen hatten Treatment-Effekt? Welche Placebo?
actice_Pos <- which(Data[,2] == 1)
placebo_Pos <- which(Data[,2] == 0)

### Wie viele jeweils?
n_active <- length(which(Data[,2] == 1))
n_placebo <- length(which(Data[,2] == 0))


### Daten extrahieren getrennt nach Active und Placebo:
Data_active <- Data[actice_Pos,]
Data_active <- as.data.frame(Data_active)
colnames(Data_active) <- varnames

Data_placebo <- Data[placebo_Pos,]
Data_placebo <- as.data.frame(Data_placebo)
colnames(Data_placebo) <- varnames


### Treatment-Effekt vorhanden?
# Grouped
bleed=c(rep(0 , 2), rep(1 , 2) , rep(2 , 2) , rep(3 , 2), rep(4 , 2),rep(5 , 2))
Treatment = rep(c("ACTIVE","PLACEBO"),6)
values_relative <- c(length(which(Data_active[,8]==0))/n_active, length(which(Data_placebo[,8]==0))/n_placebo , length(which(Data_active[,8]==1))/n_active, length(which(Data_placebo[,8]==1))/n_placebo, length(which(Data_active[,8]==2))/n_active, length(which(Data_placebo[,8]==2))/n_placebo, length(which(Data_active[,8]==3))/n_active, length(which(Data_placebo[,8]==3))/n_placebo, length(which(Data_active[,8]==4))/n_active, length(which(Data_placebo[,8]==4))/n_placebo, length(which(Data_active[,8]==5))/n_active, length(which(Data_placebo[,8]==5))/n_placebo)
data_graphics=data.frame(bleed,Treatment,values_relative)

plot1 <- ggplot(data_graphics, aes(fill=Treatment, y=values_relative, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nosebleeds (all patients)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)



### Unterschiede zwischen Gruppen aufgrund von t-Test:
t.test(Data_active[,8], Data_placebo[,8], alternative = "two.sided", mu = 0, conf.level = 0.95)


### Wirkung wird von Poissonregression bestätigt (bi- und multivariat).
summary(glm(nosebleeds ~ arm, data=Data, family = poisson))
summary(glm(nosebleeds ~ arm + country + mucus.viscosity + tissue.use + previous.year, data=Data, family = poisson))

### Overdispisperion? Negative-Binomial-Regression:
summary(glm.nb(nosebleeds ~ arm, data = Data))
summary(glm.nb(nosebleeds ~ arm + country + mucus.viscosity + tissue.use + previous.year, data = Data))




################################# How might you show how the treatment effect depends on nasal mucus viscosity?
### Zerlege Datensatz: Patienten, die zum unteren Quartil bezüglich nasal mucus viscosity gehören.
values_relative_025 <- c(length(which(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8]==0))/length(which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8]==0))/length(which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T))) , length(which(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8]==1))/length(which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8]==1))/length(which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8]==2))/length(which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8]==2))/length(which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8]==3))/length(which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8]==3))/length(which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8]==4))/length(which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8]==4))/length(which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8]==5))/length(which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8]==5))/length(which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T))))
data_graphics_025 <- data.frame(bleed,Treatment,values_relative_025)

plot2 <-ggplot(data_graphics_025, aes(fill=Treatment, y=values_relative_025, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nosebleeds (Low. 25% of nasal mucus visc.)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)


### Patienten, die zum oberen Quartil bezüglich nasal mucus viscosity gehören.
values_relative_075 <- c(length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8]==0))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8]==0))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T))) , length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8]==1))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8]==1))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8]==2))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8]==2))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8]==3))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8]==3))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8]==4))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8]==4))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8]==5))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8]==5))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T))))
data_graphics_075 <- data.frame(bleed,Treatment,values_relative_075)

plot3 <- ggplot(data_graphics_075, aes(fill=Treatment, y=values_relative_075, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nosebleeds (High. 25% of nasal mucus visc.)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)


### Patienten, die zum letztel Dezil bezüglich nasal mucus viscosity gehören.
values_relative_090 <- c(length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8]==0))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8]==0))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T))) , length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8]==1))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8]==1))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8]==2))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8]==2))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8]==3))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8]==3))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8]==4))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8]==4))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8]==5))/length(which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T))), length(which(Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8]==5))/length(which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T))))
data_graphics_090 <- data.frame(bleed,Treatment,values_relative_090)

plot4 <-  ggplot(data_graphics_090, aes(fill=Treatment, y=values_relative_090, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nosebleeds (High. 10% of nasal mucus visc.)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)

### Alle Balkendiagramme in einem Plot
grid.arrange(plot2, plot3, ncol=2)
grid.arrange(plot1, plot2, plot3, plot4, ncol=2) ### Starker Effect für höhere Perzentile!



### Unterschiede zwischen Gruppen aufgrund von t-Test:
t.test(Data_active[which(Data_active$mucus.viscosity < quantile(Data_active$mucus.viscosity, 0.25, na.rm = T)),8], Data_placebo[which(Data_placebo$mucus.viscosity < quantile(Data_placebo$mucus.viscosity, 0.25, na.rm = T)),8], alternative = "two.sided", mu = 0, conf.level = 0.95)
t.test(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.75, na.rm = T)),8], Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.75, na.rm = T)),8], alternative = "two.sided", mu = 0, conf.level = 0.95)
t.test(Data_active[which(Data_active$mucus.viscosity > quantile(Data_active$mucus.viscosity, 0.90, na.rm = T)),8], Data_placebo[which(Data_placebo$mucus.viscosity > quantile(Data_placebo$mucus.viscosity, 0.90, na.rm = T)),8], alternative = "two.sided", mu = 0, conf.level = 0.95)



### Gibt es auch den höheren Einfluss mittels Poisson-Regression?
coef(summary(glm(nosebleeds ~ arm, data=Data, family = poisson)))
coef(summary(glm(nosebleeds ~ arm, data=Data[which(Data$mucus.viscosity < quantile(Data$mucus.viscosity, 0.25, na.rm = T)),], family = poisson)))
coef(summary(glm(nosebleeds ~ arm, data=Data[which(Data$mucus.viscosity > quantile(Data$mucus.viscosity, 0.75, na.rm = T)),], family = poisson)))
coef(summary(glm(nosebleeds ~ arm, data=Data[which(Data$mucus.viscosity > quantile(Data$mucus.viscosity, 0.90, na.rm = T)),], family = poisson)))




########################## What about the effect of paper tissues? Effekt besonders hoch für Patienten mit hoher Anzahl an Tissues? Ja!!!
parts <- data.frame(
  names = c("High","Medium"),
  vals = c(length(which(Data$tissue.use == 1)), length(which(Data$tissue.use == 0)))
)

waffle(parts, rows = 15,colors = c("#fb8072", "#8dd3c7", "white"), title = "Tissue use frequencies")+theme(plot.title = element_text(hjust = 0.5))


### Unterschiede der Treatment-Effekte für Patienten mit hoem und niedrigen Verbauch?
t.test(Data_active[which(Data_active$mucus.viscosity > mean(Data_active$mucus.viscosity, na.rm = T)),8], Data_placebo[which(Data_placebo$mucus.viscosity > mean(Data_placebo$mucus.viscosity, na.rm = T)),8], alternative = "two.sided", mu = 0, conf.level = 0.95)
summary(glm(nosebleeds ~ arm, data=Data[which(Data$mucus.viscosity > mean(Data$mucus.viscosity, na.rm = T)),][which(Data$tissue.use == 1),], family = poisson))
summary(glm(nosebleeds ~ arm, data=Data[which(Data$mucus.viscosity > mean(Data$mucus.viscosity, na.rm = T)),][which(Data$tissue.use == 0),], family = poisson))



### Treatment-Effekt bei Menschen mit hohem Taschentuchverbaucht stärker als bei Medium-Verbrauch?
values_relative_high <- c(length(which(Data_active[which(Data_active$tissue.use == 1),8]==0))/nrow(Data_active[which(Data_active$tissue.use == 1),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 1),8]==0))/nrow(Data_placebo[which(Data_placebo$tissue.use == 1),]) , length(which(Data_active[which(Data_active$tissue.use == 1),8]==1))/nrow(Data_active[which(Data_active$tissue.use == 1),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 1),8]==1))/nrow(Data_placebo[which(Data_placebo$tissue.use == 1),]),length(which(Data_active[which(Data_active$tissue.use == 1),8]==2))/nrow(Data_active[which(Data_active$tissue.use == 1),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 1),8]==2))/nrow(Data_placebo[which(Data_placebo$tissue.use == 1),]), length(which(Data_active[which(Data_active$tissue.use == 1),8]==3))/nrow(Data_active[which(Data_active$tissue.use == 1),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 1),8]==3))/nrow(Data_placebo[which(Data_placebo$tissue.use == 1),]), length(which(Data_active[which(Data_active$tissue.use == 1),8]==4))/nrow(Data_active[which(Data_active$tissue.use == 1),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 1),8]==4))/nrow(Data_placebo[which(Data_placebo$tissue.use == 1),]), length(which(Data_active[which(Data_active$tissue.use == 1),8]==5))/nrow(Data_active[which(Data_active$tissue.use == 1),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 1),8]==5))/nrow(Data_placebo[which(Data_placebo$tissue.use == 1),]))
data_graphics_high=data.frame(bleed,Treatment,values_relative_high)
plot1_high <- ggplot(data_graphics_high, aes(fill=Treatment, y=values_relative_high, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nosebleeds (High tissue use)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)


### Medium-Verbrauch:
values_relative_medium <- c(length(which(Data_active[which(Data_active$tissue.use == 0),8]==0))/nrow(Data_active[which(Data_active$tissue.use == 0),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 0),8]==0))/nrow(Data_placebo[which(Data_placebo$tissue.use == 0),]) , length(which(Data_active[which(Data_active$tissue.use == 0),8]==1))/nrow(Data_active[which(Data_active$tissue.use == 0),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 0),8]==1))/nrow(Data_placebo[which(Data_placebo$tissue.use == 0),]),length(which(Data_active[which(Data_active$tissue.use == 0),8]==2))/nrow(Data_active[which(Data_active$tissue.use == 0),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 0),8]==2))/nrow(Data_placebo[which(Data_placebo$tissue.use == 0),]), length(which(Data_active[which(Data_active$tissue.use == 0),8]==3))/nrow(Data_active[which(Data_active$tissue.use == 0),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 0),8]==3))/nrow(Data_placebo[which(Data_placebo$tissue.use == 0),]), length(which(Data_active[which(Data_active$tissue.use == 0),8]==4))/nrow(Data_active[which(Data_active$tissue.use == 0),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 0),8]==4))/nrow(Data_placebo[which(Data_placebo$tissue.use == 0),]), length(which(Data_active[which(Data_active$tissue.use == 0),8]==5))/nrow(Data_active[which(Data_active$tissue.use == 0),]), length(which(Data_placebo[which(Data_placebo$tissue.use == 0),8]==5))/nrow(Data_placebo[which(Data_placebo$tissue.use == 0),]))
data_graphics_medium=data.frame(bleed,Treatment,values_relative_medium)
plot1_MEDIUM <- ggplot(data_graphics_medium, aes(fill=Treatment, y=values_relative_medium, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nodebleeds (Medium tissue use)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)

### Grafiken nebeneinander
grid.arrange(plot1_MEDIUM, plot1_high,ncol=2)



### Unterschiede in den durchschnittlichen Hospitalization Days zwischen Patienten viel und wenig Tissue Use?:
t.test(Data$nosebleeds[which(Data$tissue.use==1)],Data$nosebleeds[which(Data$tissue.use==0)])



### Ist der Effekt des Treatments in einer der Gruppen größer?
### Wie groß ist der Treatment-Effekt für Personen mit hohem Verbrauch an Taschentüchern:
t.test(Data$nosebleeds[which(Data$tissue.use==1 & Data$arm==1)],Data$nosebleeds[which(Data$tissue.use==1  & Data$arm==0)])
t.test(Data$nosebleeds[which(Data$tissue.use==0 & Data$arm==1)],Data$nosebleeds[which(Data$tissue.use==0  & Data$arm==0)])

# Bei Patienten mit hohem Verbrauch KEIN Treatment-Effekt zu beobachten!
coef(summary(glm(nosebleeds ~ arm, data=Data[which(Data$tissue.use==1),], family = poisson)))
coef(summary(glm(nosebleeds ~ arm, data=Data[which(Data$tissue.use==0),], family = poisson)))



########################### Effekt des Landes:
#################################  Hospitalization for nosebleed may depend on local medical practice. Does this have any impact? How can you understand this?

### Haben die Länder einen Effekt auf Nasenbluten-Häufigkeit?
model_nosebleed <- glm(nosebleeds ~ country , data=Data, family = poisson)

### Visualisierung:
country_levels <- levels(Data$country)
country_means <- rep(NA,length(country_levels))
for(k in 1:length(country_levels)){
  country_means[k] <- mean(Data$nosebleeds[which(Data$country == country_levels[k])])
}

country_means_former <- rep(NA,length(country_levels))
for(k in 1:length(country_levels)){
  country_means_former[k] <- mean(Data$previous.year[which(Data$country == country_levels[k])])
}
country_levels_both <- rep(country_levels ,2)
Period <- c(rep("Previous Year",10),rep("During Study",10))
country_means_both <- c(country_means_former, country_means)
plot_data_former_both <- data.frame(country_levels_both,country_means_both,Period)

### Vergleiche durchschnittliche Aufenthalte von letzten und aktuellen Jahr getrennt nach Ländern
ggplot(plot_data_former_both, aes(fill=Period, y=country_means_both, x=country_levels_both))+
  geom_point(aes(colour = Period), size = 4)+
  labs(x = "Country", y = "Arithmetic Mean")+
  ylim(0,3)+
  ggtitle("Country specific means for nosebleeds")+
  theme(plot.title = element_text(hjust = 0.5))+geom_path(aes(group = country_levels_both))



### Betrachte, ob sich die durchschnittlichen Aufenthalte nach Ländern für Active- und Placebo-Gruppe unterscheiden:
country_levels_both <- rep(country_levels ,2)
Treatment <- c(rep("ACTIVE",10),rep("PLACEBO",10))
country_means_Treatment <- rep(NA,20)
for(i in 1:10){
  country_means_Treatment[i] <- mean(Data_active[which(Data_active$country==toupper(letters[i])),8])
}
for(i in 11:20){
  country_means_Treatment[i] <- mean(Data_placebo[which(Data_placebo$country==toupper(letters[i-10])),8])
}

plot_data_former_both <- data.frame(country_levels_both,country_means_Treatment,Treatment)
ggplot(plot_data_former_both, aes(fill=Treatment, y=country_means_Treatment, x=country_levels_both))+
  geom_point(aes(colour = Treatment), size = 4)+
  labs(x = "Country", y = "Arithmetic Mean")+
  ylim(0,1.2)+
  ggtitle("Country specific means for nosebleeds")+
  theme(plot.title = element_text(hjust = 0.5))+geom_path(aes(group = country_levels_both))


### Unterschiedliche Treatment-Effekte nach Ländern?
t_tests_countries <- vector(mode="list",10)
t_tests_countries_means <- matrix(NA,10,3)
for(i in 1:10){
 t_tests_countries[[i]] <- t.test(Data_active[which(Data_active$country==toupper(letters[i])),8], Data_placebo[which(Data_placebo$country==toupper(letters[i])),8], alternative = "two.sided", mu = 0, conf.level = 0.95)
 t_tests_countries_means[i,] <- c(as.numeric(t_tests_countries[[i]]$estimate),t_tests_countries[[i]]$p.value)
}


#### Analysis of variance
variance_analysis <- summary(aov(mucus.viscosity ~ country, data = Data))


##### Ist der Effekt durch superdupripine unterschiedlich, je nachdem wie oft die Patienten
# letztes Jahr im Krankhaus waren?
t.test(Data_active[which(Data_active$previous.year>3),8], Data_placebo[which(Data_placebo$previous.year>3),8], alternative = "two.sided", mu = 0, conf.level = 0.95)
t.test(Data_active[which(Data_active$previous.year==2),8], Data_placebo[which(Data_placebo$previous.year==2),8], alternative = "two.sided", mu = 0, conf.level = 0.95)

# --> Diejenigen, die letztes Jahr "nur" zweimal im Krankenhaus waren, zeigen größeren Treatment-Effekt


################################# How might you predict the rate of nosebleed from the data that you have? What might a statistical model for this look like?
######################## Mit Zero-inflated Poisson Möglichkeit der Prediction aufzeigen:
zeroinflexample <- zeroinfl(nosebleeds ~ arm + country + mucus.viscosity + tissue.use + previous.year, data=Data)
predict(zeroinflexample, type="prob")


minCut = 5
minDev = 0.0001
treeimp <- rpart(nosebleeds ~ arm + mucus.viscosity + tissue.use, data = Data, method = "class",
                 control = rpart.control(minbucket = 5, cp = minDev))  

summary(treeimp)
printcp(treeimp)
plotcp(treeimp)
plot(treeimp, uniform=TRUE,
     main="Classification Tree for Nosebleeds")
text(treeimp, use.n=TRUE, all=TRUE, cex=.6)


################ Wie lang waren die Befragten letztes JAhr im Hospital wegen Nasenbluten`?
table(Data$previous.year)
table(Data$nosebleeds)
data_previous_now <- c(length(which(Data$previous.year==0))/n,length(which(Data$nosebleeds==0))/n,length(which(Data$previous.year==1))/n,length(which(Data$nosebleeds==1))/n,length(which(Data$previous.year==2))/n,length(which(Data$nosebleeds==2))/n,length(which(Data$previous.year==3))/n,length(which(Data$nosebleeds==3))/n,length(which(Data$previous.year==4))/n,length(which(Data$nosebleeds==4))/n,length(which(Data$previous.year==5))/n,length(which(Data$nosebleeds==5))/n)
Year <- rep(c("Previous Year","During study"),6)
bleed=c(rep(0 , 2), rep(1 , 2) , rep(2 , 2) , rep(3 , 2), rep(4 , 2),rep(5 , 2)) 
data_graphics_0 <- data.frame(bleed,Year,data_previous_now)


ggplot(data_graphics_0, aes(fill=Year, y=data_previous_now, x=bleed)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Amount of Nosebleeds", y = "Relative Frequency") +
  ggtitle("Nosebleeds during study and previous year") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = 0:5)+ 
  ylim(0,1)

# --> Dieses Jahr (während der Behandlung) viel seltener.