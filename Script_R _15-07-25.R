##The effects of agricultural landscape on the functional diversity of bird species in the Neotropical# #dryland                                                                                              #
#Mônica da Costa Lima, Thiago Nascimento Zanetti, Helder Farias Pereira de Araujo                     #
#                                                                                                     #
## SCRIPT FUNCTIONAL DIVERSITY (FD), GLMS AND PCoA                                                    #     
#Diversidade funcional no R baseados em scripts de Julia C. Sfair                                     #
#######################################################################################################

#a ordem das especies na matriz de comunidade deve ser EXATAMENTE
#a mesma da ordem das especies da matriz de atributos. 

##DEFINIR DIRETORIO DE TRABALHO
setwd ("WORK DIRECTORY")
dir()

#A matriz deve ser assim:
#- Matriz de comunidade:
#linhas sao areas e colunas sao especies-
#Matriz de atributos: linhas sao especies e colunas sao atributos 

#ABRIR PLANILHAS
comm <- read.csv("comm_abundancias.txt", sep = "\t", header = T, row.names = 1)

att <- read.table ("att para 90 sps.txt", sep = "\t", h = T, row.names = 1)

##TRANSFORMAR TODAS AS COLUNAS EM OBJETOS
attach (att)

colnames(comm) <- rownames (att)
att

##INSTALAR- SE NECESSARIO- E ABRIR OS PACOTES
install.packages("FD")
library(FD)
library (vegan)

#Para se calcular o indice FD, de Petchey & Gaston (2002. Functional 
#diversity (FD), species richness and community composition. Ecol. Lett. 
#5, 402 411.

library (dplyr)
att %>%
  dplyr::summarise_all(class) %>%
  tidyr::gather(variable, class)

#######################
#FD
######################
## 3. Combinar cada conjunto de atributos de acordo com sua natureza em um
# data.frame separado.
# 3.1. interger- ordinais # foram reconhecidas como numéricas ou categóricas.
names (att)

att$terrestre<- as.ordered(att$Terrestre)
att$Sub.bosque <- as.ordered(att$Sub.bosque)
att$Médio_Sub.bosque<- as.ordered(att$Médio_Sub.bosque)
att$dossel <- as.ordered(att$Dossel)
att$Areial <- as.ordered(att$Areial)
att$agua <- as.ordered(att$Água)
att$vertebrados <- as.ordered(att$Vertebrados)
att$invertebrados <- as.ordered(att$Invertebrados)
att$Material_vegetal <- as.ordered(att$Material_vegetal)
att$Frutos <- as.ordered(att$Frutos)
att$Grãos <- as.ordered(att$Grãos)
att$Nectar <- as.ordered(att$Nectar)
att$Saprofago <- as.ordered(att$Saprofago)

att_ord <- cbind.data.frame(terrestre = att$Terrestre,
                            Sub.bosque = att$Sub.bosque,
                            Médio_Sub.bosque = att$Médio_Sub.bosque, dossel = att$dossel, Areial = att$Areial, agua= att$Água, vertebrados= att$Vertebrados, invertebrados= att$Invertebrados, Material_vegetal= att$Material_vegetal, Frutos= att$Frutos, Grãos= att$Grãos, Nectar= att$Nectar, Saprofago= att$Saprofago)

rownames(att_ord) <- rownames(att)

# Agora, combinar os dois data.frames em uma lista chamada "ktab".
massa= data.frame(att$massa.g.)#precisa ser data frame
ktab_list <- ktab.list.df(list(att_ord, massa))

# Por fim, calcular a distância funcional entre as espécies.
# Em "type", a letra "N" indica variável categórica (ou nominal),
# enquanto a letra "O" indica variável ordinal.
dist_mist <- dist.ktab(ktab_list, type = c("O", "Q"))

## Visualize os dados com uma PCoA
library(ggplot2)
library(ggrepel)

#normalizar
colnames(comm) <- rownames (att.norm)
att.norm <- decostand(dist_mist, "normalize", MARGIN = 2,  na.rm = T)

nomes_comm <- colnames(comm)
nomes_dist_mist <- rownames(dist_mist)
diferenca <- setdiff(nomes_comm, nomes_dist_mist)

dist.att <- gowdis(att.norm)
agrup <- (hclust(dist.att , "average")) #UPGMA
plot(hclust(dist.att, "average"))
treedive(comm, tree = agrup, match.force = TRUE)

FD <- treedive(comm, tree = agrup, match.force = TRUE)
treeheight(agrup)

## SALVAR MATRIZ
write.csv(FD, "FD_results.gower.csv")

##OUTROS INDICES DE DIVERSIDADE FUNCIONAL
# Diversidade funcional
install.packages("ade4")
install.packages("vegan")
install.packages("lawstat") # pacote de normalidade
library(ade4)
library(vegan)
library(lawstat) 
#install.packages("FD")
library(FD)

# Agora, combinar os dois data.frames em uma lista chamada "ktab".

#att1<-read.table("att para 90 sps.txt",sep = "\t", header = T, row.names = 1)

att1 <- read.table("att.txt", sep = "\t", header = T,row.names = 1)
str (dist_mist)
comm <- read.csv("comm_pool.csv", sep = ";", header = T, row.names = 1)
colnames(comm) <- rownames (att1)

#Composição funcional
functcomp(x = att, a = as.matrix(comm), CWM.type = c("dom"))
functcomp(x = att, a = as.matrix(comm), CWM.type = c("all"))

#Normalização
att.norm <- decostand(att, "normalize", MARGIN = 2,  na.rm = T)

dist.att <- gowdis(att.norm)
plot(hclust(dist.att, "average"))
agrup <- (hclust(dist.att, "average"))
#Testando correlação cofenética
cophenetic(agrup)
cor(dist.att, cophenetic(agrup))
#PCA para verificar correlação entre variáveis
att1= att[,-1]
biplot(prcomp(att1, scale = T))
#Para ver os escores das espécies
prcomp(att1, scale = T)
#Para ver a porcentagem de explicação de cada eixo
summary(prcomp(att1, scale = T))
#diversidade funcional
div.func <- dbFD(x = dist_mist, a = comm, w.abun = T, stand.x = T, stand.FRic = T, m = 5)#Para rodar essa função com a matriz de distancia ktab é so transformar ela em matriz com a função as.matrix.
ncol (att1)
#É importante ler o help da função acima, pois ela calcula vários índices e costuma ter problemas.

#Para ver os resultados
div.func

#Para ver os índices em separado
div.func$FRic # riqueza funcional
div.func$FEve # uniformidade funcional

write.csv (div.func$FRic, "div.func$FRic.csv")
write.csv (div.func$FEve, "div.func$FEve.csv")
# dispersão funcional, 
fdisp(d = gowdis(att), a = as.matrix(comm))

div.func$CWM # dominancia de traits
div.func$RaoQ # entropia de Rao
##############################

##GLMS
##############################
##PARA OS MODELOS E PRECISO ABRIR A PLANILHA COM OS DADOS DE PAISAGEM E DE FD
##############################
#corrplot
###############################
setwd ("WORK DIRECTORY")
getwd ()
dir ()
corr <- read.table("Het-3.txt", sep = "\t", header = T, row.names = 1)
attach (corr)

library (corrplot)

#tiff(filename = "Corr.FD-landscape.tiff",width = 800, height = 600, units = "px", pointsize = 20)
corrplot(corr=cor(corr[,]),type="upper", order="hclust",method = "number")

dev.copy2pdf (width=15,height=15, file="corrplot 23-04-24.pdf")


###############################
##ESCOLHER DIRETORIO DE TRABALHO
setwd ("WORK DIRECTORY")
getwd ()
dir ()
heterogeneidade <- read.table("Het-3.txt", sep = "\t", header = T, row.names = 1)
attach (heterogeneidade)
library (ggplot2)

tiff(filename = "Comple-FD.tiff",width = 1000, height = 800, units = "px", pointsize = 20)
ggplot(data=heterogeneidade,aes(Compl.nat, FD
))+geom_point()+ theme_classic ()+  labs( x = "Landscape Complexity",
                                          y = "FD") +geom_smooth(method="lm")
dev.off()

tiff(filename = "Compo-FD.fiff",width = 1000, height = 800, units = "px", pointsize = 20)
ggplot(data=heterogeneidade,aes(Compo.nat, FD
))+geom_point()+ theme_classic ()+  labs( x = "Landscape Composition",
                                          y = "FD") +geom_smooth(method="lm")
dev.off()

tiff(filename = "Comple-FRic.fiff",width = 1000, height = 800, units = "px", pointsize = 20)
ggplot(data=heterogeneidade,aes(Compl.plot, FRic
))+geom_point()+ theme_classic ()+  labs( x = "Landscape Complexity",
                                          y = "FRic") +geom_smooth(method="lm")
dev.off()

tiff(filename = "Compo-FRic.fiff",width = 1000, height = 800, units = "px", pointsize = 20)
ggplot(data=heterogeneidade,aes(Compo.nat, FRic
))+geom_point()+ theme_classic ()+  labs( x = "Landscape Composition",
                                          y = "FRic") +geom_smooth(method="lm")
dev.off()

#####
#figura completa
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(gridExtra)

# Função para formatar os rótulos dos eixos
formato_eixo <- function(grafico, remove_x = FALSE, remove_y = FALSE) {
  grafico + 
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 16), # Tamanho do título do eixo X
      axis.title.y = element_text(size = 16), # Tamanho do título do eixo Y
      axis.text.x = if(remove_x) element_blank() else element_text(size = 12),  # Remover ou ajustar rótulos do eixo X
      axis.text.y = if(remove_y) element_blank() else element_text(size = 12)   # Remover ou ajustar rótulos do eixo Y
    )
}

# Criar as figuras individuais
fig1 <- ggplot(data = heterogeneidade, aes(Compl.plot, FD)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "", y = "FD")

fig2 <- ggplot(data = heterogeneidade, aes(Compo.nat, FD)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "", y = "")

# Aplicar formato aos rótulos dos eixos com remoção específica
fig1 <- formato_eixo(fig1, remove_x = TRUE)
fig2 <- formato_eixo(fig2, remove_x = TRUE, remove_y = TRUE)

formato_eixo <- function(grafico, remove_x = FALSE, remove_y = FALSE) {
  grafico + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 16), # Tamanho do título do eixo X
      axis.title.y = element_text(size = 16), # Tamanho do título do eixo Y
      axis.text.x = if(remove_x) element_blank() else element_text(size = 12),  # Remover ou ajustar rótulos do eixo X
      axis.text.y = if(remove_y) element_blank() else element_text(size = 12)   # Remover ou ajustar rótulos do eixo Y
    )
}

fig3 <- ggplot(data = heterogeneidade, aes(Compl.plot, FRic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Landscape Complexity", y = "FRic")

fig4 <- ggplot(data = heterogeneidade, aes(Conf.nat, FRic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Landscape Configuration", y = "FEve")

# Aplicar formato aos rótulos dos eixos com remoção específica
fig3 <- formato_eixo(fig3)
fig4 <- formato_eixo(fig4)

# Combinar os gráficos em uma única figura 2 por 2
combined_plots <- grid.arrange(fig1, fig2, fig3, fig4, ncol = 2)

#tiff(filename = "Fig 2 regressoes.tiff",width = 1200, height = 800, units = "px", pointsize = 300)
#combined_plots <- grid.arrange(fig1, fig2, fig3, fig4, ncol = 2)
#dev.off()

#salvar em pdf
dev.copy2pdf (width=8,height=5, file="fig 2 regre 24-06.pdf")

#tiff(filename = "ConfNat-FEve.fiff",width = 1000, height = 800, units = "px", pointsize = 20)
#ggplot(data=heterogeneidade,aes(Conf.nat, FEve
#))+geom_point()+ theme_classic ()+  labs( x = "Natural cover configuration",
 #                                         y = "FEve") +geom_smooth#(method="lm")
#dev.off()

#### plots 2024
# Library
library(ggplot2)
library(hrbrthemes)

# linear trend + confidence interval
p3 <- ggplot(heterogeneidade, aes(x=Compl.plot, y=my_y)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

##################
##GLMs
####################
## PACOTES NECESSARIOS
library(nlme)
library (MuMIn)

##MODELOS COM FD RETIRANDO VARIAVIES
#Compl.plot+  Compo.nat+ Conf.nat

GLM_FD<-glm(FD~Compl.plot+  Conf.nat, data = heterogeneidade)
summary(GLM_FD)
#AICc(Full_GLM, k = 2, REML = NULL)

backwards = step(GLM_FD, direction="both")

#####
GLM_FRic<-glm(FRic~Compl.plot, data = heterogeneidade)
summary(GLM_FRic)
#AICc(GLM_2, k = 2, REML = NULL)

backwards = step(GLM_FRic, direction="both")
#####

GLM_FEve<-glm(FEve~  Conf.nat, data = heterogeneidade)
summary(GLM_FEve)
#AICc(GLM_FEve, k = 2, REML = NULL)

backwards = step(GLM_FEve, direction="both")

## AQUI TESTEI AS DIFERENTES VARIAVEIS JUNTO A COMPLEXIDADE

##GLM PARA CADA VARIAVEL SEPARADAMENTE
#####
##Compo.nat
GLM_COMPO_NAT<-glm(FD~Compo.nat, data = heterogeneidade)
summary(GLM_COMPO_NAT)
AICc(GLM_COMPO_NAT, k = 2, REML = NULL)

##Conf.nat
GLM_CONF_NAT<-glm(FD~Conf.nat, data = heterogeneidade)
summary(GLM_CONF_NAT)
AICc(GLM_CONF_NAT, k = 2, REML = NULL)

##PLOTS
par(mfrow = c(2,3))

#COMPLEXIDADE
plot(FRic~Compl.plot, ylab = "", xlab="Landscape complexity", col="black", data = heterogeneidade)
abline(lm(FRic~Compl.plot),col="black")

#COMPOSI  O NATURAL
plot(FRic~Compo.nat, ylab = "", xlab="Natural cover composition", col="black", data = heterogeneidade)
abline(lm(FRic~Compo.nat),col="black")

#COMPOSI  O IMPACTADA
plot(FRic~ Compo.impac, ylab = "", xlab="Shrubland cover composition", col="black", data = heterogeneidade)
abline(lm(FRic~Compo.impac),col="black")

#CONFIGURA  O NATURAL
plot(FRic~Conf.nat, ylab = "", xlab="Natural cover configuration", col="black", data = heterogeneidade)
abline(lm(FRic~Conf.nat),col="black")

#CONFIGURA  O IMPACTADA
plot(FRic~conf.impac, ylab = "", xlab=" Shrubland cover configuration", col="black", data = heterogeneidade)
abline(lm(FRic~conf.impac),col="black")

#CONFIGURA  O CROP
plot(FRic~conf.crop, ylab = "", xlab="Crop cover configuration", col="black", data = heterogeneidade)
abline(lm(FRic~conf.crop),col="black")

##REGRESS O LINEAR
plot(FRic~Compl.plot, ylab = "Functional richness", xlab="Landscape complexity", col="black",
data = heterogeneidade, las=1, family="sans", pch=16, cex.lab=1.2, cex.axis=0.9)
text(8.8,2.2,"R =0.46", cex=0.8)
text(8.8,2,"p<0.01", cex=0.8)
abline(lm(FEve~Compl.plot),col="black")
#################################
##GLM com parametros solicitados pela revista
# Ajustar o modelo GLM
################
Full_GLM <- glm(FD ~ Compl.plot + Compo.nat+ Conf.nat , data = heterogeneidade)

# Resumo do modelo GLM
summary(Full_GLM)

# Obtendo o número de parâmetros
num_parametros <- length(coef(Full_GLM))

# Imprimindo o resultado
print(num_parametros)

# Calcular AICc (Akaike Information Criterion com correção para amostras pequenas)
library(AICcmodavg)
aic <- AICc(Full_GLM, k = 2, REML = FALSE)

# Calcular delta AICc e AICc weight
delta_aicc <- AICc(Full_GLM, k = 2, REML = FALSE) - min(aic)
aic_weight <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Imprimir AICc, delta AICc e AICc weight
cat("AICc:", aic, "\n")
cat("Delta AICc:", delta_aicc, "\n")
cat("AICc Weight:", aic_weight, "\n")

# Log-likelihood
log_likelihood <- logLik(Full_GLM)

# Extraindo estimativas e precisão (erros padrão e intervalo de confiança)
estimates <- coef(Full_GLM)
standard_errors <- summary(Full_GLM)$coefficients[, "Std. Error"]
conf_interval <- confint(Full_GLM)

# Imprimir log-likelihood, estimativas e precisão
cat("Log-Likelihood:", log_likelihood, "\n")
cat("Parameter Estimates:\n")
print(estimates)
cat("Standard Errors:\n")
print(standard_errors)
cat("Confidence Intervals:\n")
print(conf_interval)

##Complexity model
GLM <- glm(FD ~ Compl.plot, data = heterogeneidade)

# Resumo do modelo GLM
summary(GLM)

# Obtendo o número de parâmetros
num_parametros <- length(coef(GLM))

# Imprimindo o resultado
print(num_parametros)

# Calcular AICc (Akaike Information Criterion com correção para amostras pequenas)
library(AICcmodavg)
aic <- AICc(GLM, k = 2, REML = FALSE)

# Calcular delta AICc e AICc weight
delta_aicc <- AICc(GLM, k = 2, REML = FALSE) - min(aic)
aic_weight <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Imprimir AICc, delta AICc e AICc weight
cat("AICc:", aic, "\n")
cat("Delta AICc:", delta_aicc, "\n")
cat("AICc Weight:", aic_weight, "\n")

# Log-likelihood
log_likelihood <- logLik(GLM)

# Extraindo estimativas e precisão (erros padrão e intervalo de confiança)
estimates <- coef(GLM)
standard_errors <- summary(GLM)$coefficients[, "Std. Error"]
conf_interval <- confint(GLM)

# Imprimir log-likelihood, estimativas e precisão
cat("Log-Likelihood:", log_likelihood, "\n")
cat("Parameter Estimates:\n")
print(estimates)
cat("Standard Errors:\n")
print(standard_errors)
cat("Confidence Intervals:\n")
print(conf_interval)

##Full GLM FRic
Full_GLM <- glm(FRic ~ Compl.plot +Compo.nat +Conf.nat, data = heterogeneidade)

# Resumo do modelo GLM
summary(Full_GLM)

# Obtendo o número de parâmetros
num_parametros <- length(coef(Full_GLM))

# Imprimindo o resultado
print(num_parametros)

# Calcular AICc (Akaike Information Criterion com correção para amostras pequenas)
library(AICcmodavg)
aic <- AICc(Full_GLM, k = 2, REML = FALSE)

# Calcular delta AICc e AICc weight
delta_aicc <- AICc(Full_GLM, k = 2, REML = FALSE) - min(aic)
aic_weight <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Imprimir AICc, delta AICc e AICc weight
cat("AICc:", aic, "\n")
cat("Delta AICc:", delta_aicc, "\n")
cat("AICc Weight:", aic_weight, "\n")

# Log-likelihood
log_likelihood <- logLik(Full_GLM)

# Extraindo estimativas e precisão (erros padrão e intervalo de confiança)
estimates <- coef(Full_GLM)
standard_errors <- summary(Full_GLM)$coefficients[, "Std. Error"]
conf_interval <- confint(Full_GLM)

# Imprimir log-likelihood, estimativas e precisão
cat("Log-Likelihood:", log_likelihood, "\n")
cat("Parameter Estimates:\n")
print(estimates)
cat("Standard Errors:\n")
print(standard_errors)
cat("Confidence Intervals:\n")
print(conf_interval)

##Complexity model FRic

GLM1 <- glm(FRic ~ Compl.plot, data = heterogeneidade)

# Resumo do modelo GLM
summary(GLM1)

# Obtendo o número de parâmetros
num_parametros <- length(coef(GLM1))

# Imprimindo o resultado
print(num_parametros)

# Calcular AICc (Akaike Information Criterion com correção para amostras pequenas)
library(AICcmodavg)
aic <- AICc(GLM1, k = 2, REML = FALSE)

# Calcular delta AICc e AICc weight
delta_aicc <- AICc(GLM1, k = 2, REML = FALSE) - min(aic)
aic_weight <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Imprimir AICc, delta AICc e AICc weight
cat("AICc:", aic, "\n")
cat("Delta AICc:", delta_aicc, "\n")
cat("AICc Weight:", aic_weight, "\n")

# Log-likelihood
log_likelihood <- logLik(GLM1)

# Extraindo estimativas e precisão (erros padrão e intervalo de confiança)
estimates <- coef(GLM1)
standard_errors <- summary(GLM1)$coefficients[, "Std. Error"]
conf_interval <- confint(GLM1)

# Imprimir log-likelihood, estimativas e precisão
cat("Log-Likelihood:", log_likelihood, "\n")
cat("Parameter Estimates:\n")
print(estimates)
cat("Standard Errors:\n")
print(standard_errors)
cat("Confidence Intervals:\n")
print(conf_interval)

##FEve Conf.nat
Full_GLM <- glm(FEve ~ Compl.plot +Compo.nat+Conf.nat, data = heterogeneidade)

# Resumo do modelo GLM
summary(Full_GLM)

# Obtendo o número de parâmetros
num_parametros <- length(coef(Full_GLM))

# Imprimindo o resultado
print(num_parametros)

# Calcular AICc (Akaike Information Criterion com correção para amostras pequenas)
library(AICcmodavg)
aic <- AICc(Full_GLM, k = 2, REML = FALSE)

# Calcular delta AICc e AICc weight
delta_aicc <- AICc(Full_GLM, k = 2, REML = FALSE) - min(aic)
aic_weight <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Imprimir AICc, delta AICc e AICc weight
cat("AICc:", aic, "\n")
cat("Delta AICc:", delta_aicc, "\n")
cat("AICc Weight:", aic_weight, "\n")

# Log-likelihood
log_likelihood <- logLik(Full_GLM)

# Extraindo estimativas e precisão (erros padrão e intervalo de confiança)
estimates <- coef(Full_GLM)
standard_errors <- summary(Full_GLM)$coefficients[, "Std. Error"]
conf_interval <- confint(Full_GLM)

# Imprimir log-likelihood, estimativas e precisão
cat("Log-Likelihood:", log_likelihood, "\n")
cat("Parameter Estimates:\n")
print(estimates)
cat("Standard Errors:\n")
print(standard_errors)
cat("Confidence Intervals:\n")
print(conf_interval)
####################

GLM1 <- glm(FEve ~ Conf.nat, data = heterogeneidade)

# Resumo do modelo GLM
summary(GLM1)

# Calcular AICc (Akaike Information Criterion com correção para amostras pequenas)
library(AICcmodavg)
aic <- AICc(GLM1, k = 2, REML = FALSE)

# Calcular delta AICc e AICc weight
delta_aicc <- AICc(GLM1, k = 2, REML = FALSE) - min(aic)
aic_weight <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))

# Imprimir AICc, delta AICc e AICc weight
cat("AICc:", aic, "\n")
cat("Delta AICc:", delta_aicc, "\n")
cat("AICc Weight:", aic_weight, "\n")

# Log-likelihood
log_likelihood <- logLik(GLM1)

# Extraindo estimativas e precisão (erros padrão e intervalo de confiança)
estimates <- coef(GLM1)
standard_errors <- summary(GLM1)$coefficients[, "Std. Error"]
conf_interval <- confint(GLM1)

# Imprimir log-likelihood, estimativas e precisão
cat("Log-Likelihood:", log_likelihood, "\n")
cat("Parameter Estimates:\n")
print(estimates)
cat("Standard Errors:\n")
print(standard_errors)
cat("Confidence Intervals:\n")
print(conf_interval)

##Script PCoA##
################################################
#Script baseado em <https://canal6.com.br/livros_loja/Ebook_Analises_Ecologicas_no_R.pdf>

## Pacotes
library(ade4)
#library(ecodados)
library(tidyverse)
library(vegan)
library(pvclust)
#library(BiodiversityR)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(ape)
library(FactoMineR)
library(factoextra)
library(FD)
library(palmerpenguins)
library(GGally)
library(fields)
library(ade4)
#library(ggord)
library(udunits2)
library(adespatial)
library(spdep)
library(mvabund)
library(reshape)

setwd ("WORK DIRECTORY")
dir()

#ABRIR PLANILHA
Dados <- read.table("ponto por atribultos funcionais agrupados.txt", sep = "\t", header = T)
Dados

#Vamos primeiramente padronizar dos dados, calcular uma matriz de distância com método BrayCurtis e depois calcular a PCoA
## Padronização dos dados com Hellinger
dados= Dados [-1]
Dados.hel <- decostand(x = dados, method = "hellinger") 

## Cálculo da matriz de distância com método Bray-Curtis
sps.dis <- vegdist(x = Dados.hel, method = "bray") 

## PCoA
pcoa.sps <- pcoa(D = sps.dis, correction = "cailliez")

## Porcentagem de explicação do Eixo 1
explicacao_eixo1= 100 * (pcoa.sps$values[, 1]/pcoa.sps$trace)[1]

## Porcentagem de explicação dos Eixo 2
explicacao_eixo2= 100 * (pcoa.sps$values[, 1]/pcoa.sps$trace)[2]

## Porcentagem de explicação acumulada dos dois primeiros eixos
sum(100 * (pcoa.sps$values[, 1]/pcoa.sps$trace)[1:2])

## Selecionar os dois primeiros eixos
eixos <- pcoa.sps$vectors[, 1:2]

## Juntar com algum dado categórico de interesse para fazer a figura
pcoa.dat <- data.frame(pontos = Dados$Ponto, eixos)

a= read.table ("clipboard", h=T)
pcoa.dat <- data.frame(Landscape = a , eixos)

### Gráfico biplot da PCoA
ggplot(pcoa.dat, aes(x = Axis.1, y = Axis.2, fill = Landscape,
                     color = Landscape, shape = Landscape)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("darkorange", "purple4","darkgreen")) +
  labs(x = "PCO 1 (36.23%)", y = "PCO 2 (34.24%)") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic()

arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="black", length=0.15, lwd=1.5)

## Sem formas diferentes
ggplot(pcoa.dat, aes(x = Axis.1, y = Axis.2, fill = Landscape, color = Landscape)) +
  geom_point(size = 2, alpha = 0.7)  +
  scale_color_manual(values = c("darkorange", "purple4", "darkgreen")) +
  scale_fill_manual(values = c("darkorange", "purple4","darkgreen")) +
  labs(x = "PCO 1 (36.23%)", y = "PCO 2 (34.24%)") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic()


###setas
# Adicionando setas
# Plot PCoA com ggplot
library(ggplot2)
library(ggrepel)

ggplot(pcoa.dat, aes(x = Axis.1, y = Axis.2, fill = Landscape, color = Landscape)) +
  geom_point(size = 2, alpha = 0.7)  +
  scale_color_manual(values = c("darkorange", "purple4", "darkgreen")) +
  scale_fill_manual(values = c("darkorange", "purple4","darkgreen")) +
  labs(x = "PCO 1 (36.23%)", y = "PCO 2 (34.24%)") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  
  #salvar em pdf
  dev.copy2pdf (width=8,height=6, file="PCoA_Mônica_Lima-27-9-23.pdf")
  
  biplot(pcoa.sps, Dados.hel, main="",xlab=paste("PCoA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCoA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=1.5, las=1, asp=1,bty="L",lwd=0.5, family="sans")

##INSERIR LINHAS CINZA DO CENTRO
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")

dev.copy2pdf (width=10,height=8, file="setas.pdf")

#Figura composta a partir das duas acima.

## Pacotes
library(ade4)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library (ape)

# Defina o diretório de trabalho
setwd("WORK DIRECTORY")
dir()

# Abra a planilha
Dados <- read.table("ponto por atribultos funcionais agrupados.txt", sep = "\t", header = TRUE)
dados <- Dados[-1]

# Padronização dos dados com Hellinger
Dados.hel <- decostand(x = dados, method = "hellinger")

# Cálculo da matriz de distância com método Bray-Curtis
sps.dis <- vegdist(x = Dados.hel, method = "bray")

# PCoA
pcoa.sps <- pcoa(D = sps.dis, correction = "cailliez")

# Porcentagem de explicação dos Eixos
explicacao_eixo1 <- 100 * (pcoa.sps$values[, 1] / pcoa.sps$trace)[1]
explicacao_eixo2 <- 100 * (pcoa.sps$values[, 1] / pcoa.sps$trace)[2]

# Selecionar os dois primeiros eixos
eixos <- pcoa.sps$vectors[, 1:2]
colnames(eixos) <- c("Axis.1", "Axis.2")

# Juntar com algum dado categórico de interesse para fazer a figura
a <- read.table("clipboard", header = TRUE) # coluna dos grupos, os tipos de paisagem
if (nrow(a) != nrow(eixos)) {
  stop("O número de linhas em 'a' não corresponde ao número de linhas em 'eixos'")
}

pcoa.dat <- data.frame(Landscape = a$Landscape, eixos) # Assume que 'Landscape' é a coluna de interesse

# Cálculo das setas
fit <- envfit(pcoa.sps$vectors, dados)
arrows_data <- as.data.frame(fit$vectors$arrows)
arrows_data$labels <- rownames(arrows_data)
colnames(arrows_data) <- c("Axis.1", "Axis.2", "labels")

# Fator de escala para ajustar a proporção das setas
scale_factor <- 0.5 # Ajuste este valor conforme necessário para reduzir ou aumentar o tamanho das setas

# Aplicar o fator de escala às coordenadas das setas
arrows_data$Axis.1 <- arrows_data$Axis.1 * scale_factor
arrows_data$Axis.2 <- arrows_data$Axis.2 * scale_factor

# Gráfico biplot da PCoA
#tiff(filename = "Fig 3 pcoa.tiff",width = 1300, height = 700, units = "px", pointsize = 300)
ggplot(pcoa.dat, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(fill = Landscape, color = Landscape, shape = Landscape), size = 4, alpha = 0.7) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("darkorange", "purple4", "darkgreen")) +
  labs(x = paste("PCO 1 (", round(explicacao_eixo1, 2), "%)", sep = ""),
       y = paste("PCO 2 (", round(explicacao_eixo2, 2), "%)", sep = "")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  geom_segment(data = arrows_data, aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2), 
               arrow = arrow(length = unit(0.1, "cm")), color = "black") +
  geom_text_repel(data = arrows_data, aes(x = Axis.1, y = Axis.2, label = labels), 
                  color = "black", size = 5) +
  theme(axis.title.x = element_text(size = 14),   # Tamanho dos rótulos dos eixos
        axis.title.y = element_text(size = 14),   # Tamanho dos rótulos dos eixos
        axis.text.x = element_text(size = 12),    # Tamanho dos valores dos eixos
        axis.text.y = element_text(size = 12),    # Tamanho dos valores dos eixos
        legend.text = element_text(size = 12),    # Tamanho do texto da legenda
        legend.title = element_text(size = 14))   # Tamanho do título da legenda




#dev.off()

# Salvar o gráfico
ggsave("pcoa setas artigo.pdf", width = 10, height = 8)

