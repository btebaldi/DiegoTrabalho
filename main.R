
# Setup -------------------------------------------------------------------
rm(list = ls())

library(readxl)
library(ggplot2)
library(dplyr)
library(magrittr)
library(lubridate)

#' Função que testa a raiz unitaria para variaveis, de uma base de dados
#'
#' @param dados Base de dados a ser analizada
#' @param colunas Colunas a serem analizadas 
#' @param critical valor critico da analise
#' @param InfoCriteria Criterio de selecao das variaveis 
#' @param N.lags Numero de lags do ADF
#'
#' @return Tabela informando se existe ou nao raiz unitaria
#' @export 
#'
#' @examples Testa.RaizUnitaria(data, c("col1", "col2"))
#'
Testa.RaizUnitaria <- function(dados,
                               colunas=NULL,
                               critical = c("5pct", "1pct", "10pct"),
                               InfoCriteria = c("Fixed", "BIC", "AIC"),
                               N.lags = 1){
  
  # Load arguments
  InfoCriteria <- "Fixed"
  p.critical <-  c("5pct" = 0.05,  "1pct"=0.01, "10pct"=0.10)
  if(is.null(colunas)){
    colunas <- colnames(dados)
  }
  
  critical <- match.arg(critical)
  
  # biblioteca requerida
  require(urca)
  
  dados.ret <- data.frame(variavel=colunas,
                          Unit.Root = as.logical(NA),
                          Criterio = as.numeric(NA),
                          trend.tau3 = as.numeric(NA),
                          trend.phi3 = as.numeric(NA),
                          trend.phi2 = as.numeric(NA),
                          drift.tau2 = as.numeric(NA),
                          drift.phi1 = as.numeric(NA),
                          none.tau1 = as.numeric(NA),
                          
                          critical.trend.tau3 = as.numeric(NA),
                          critical.trend.phi3 = as.numeric(NA),
                          critical.trend.phi2 = as.numeric(NA),
                          critical.drift.tau2 = as.numeric(NA),
                          critical.drift.phi1 = as.numeric(NA),
                          critical.none.tau1 = as.numeric(NA)
  )
  
  for(i in 1:nrow(dados.ret)){
    # Retira os dados da tabela
    x <- dados[[dados.ret$variavel[i]]]
    
    # se existe NA na serie, considero como erro zerado.
    if(sum(is.na(x)) > 0){
      cat(sprintf("\tSerie %d contem NA\n", i));
      stop("Encontrado NA")
    }
    
    if(!is.numeric(x)) {
      cat(sprintf("\tSerie %s nao e numerica. Skip!\n", colunas[i]));
      next()
    }
    
    if(sum(abs(x - mean(x))) == 0) {
      cat(sprintf("\tSerie %s nao tem variacao. Skip!\n", colunas[i]));
      next()
    }
    
    # Estima os modelos com tendencia, drift, e none
    ADF.none <- ur.df(x, type = "none", lags = N.lags, selectlags = InfoCriteria)
    ADF.drift <- ur.df(x, type = "drift", lags = N.lags, selectlags = InfoCriteria)
    ADF.trend <- ur.df(x, type = "trend", lags = N.lags, selectlags = InfoCriteria)
    
    # Guarda as estatisiticas
    dados.ret[i, c("trend.tau3", "trend.phi3", "trend.phi2")] <- ADF.trend@teststat["statistic", c("tau3", "phi3", "phi2")]
    dados.ret[i, c("drift.tau2", "drift.phi1")] <- ADF.drift@teststat["statistic", c("tau2", "phi1")]
    dados.ret[i, c("none.tau1")] <- ADF.none@teststat["statistic", c("tau1")]
    
    
    # Guarda os Valores criticos estatisiticas
    dados.ret[i, c("critical.trend.tau3", "critical.trend.phi3", "critical.trend.phi2")] <- ADF.trend@cval[,critical]
    dados.ret[i, c("critical.drift.tau2", "critical.drift.phi1")] <- ADF.drift@cval[,critical]
    dados.ret[i, c("critical.none.tau1")] <- ADF.none@cval[,critical]
    
    
    
    if(ADF.trend@teststat["statistic", "tau3"] < ADF.trend@cval["tau3", critical]){
      # se temos rejeicao da nula entao nao tem os unit root
      dados.ret[i, "Unit.Root"] <- FALSE
      dados.ret[i, "Criterio"] <- 1
    } else {
      if(ADF.trend@teststat["statistic", "phi3"] < ADF.trend@cval["phi3", critical]) {
        # aqui a trend nao ? significante, entao devemos estimar o modelo drift
        
        if(ADF.drift@teststat["statistic", "tau2"] < ADF.drift@cval["tau2", critical]){
          #  temos rejeicao da nula e com isso nao temos raiz unit?ria
          dados.ret[i, "Unit.Root"] <- FALSE
          dados.ret[i, "Criterio"] <- 2
        } else {
          # temos de testar a presenca de drift nao significante
          if(ADF.drift@teststat["statistic", "phi1"] < ADF.drift@cval["phi1", critical]){
            # aqui nao temos drift, temos de estimar o modelo sem nada
            if(ADF.none@teststat["statistic", "tau1"] < ADF.none@cval["tau1", critical]){
              #  temos rejeicao da nula e com isso nao temos raiz unitaria
              dados.ret[i, "Unit.Root"] <- FALSE
              dados.ret[i, "Criterio"] <- 3
            } else {
              dados.ret[i, "Unit.Root"] <- TRUE
              dados.ret[i, "Criterio"] <- 4
            }
          } else {
            # temos presenca de drift. Realizamos o teste utilizando a distribuicao normal.
            if(ADF.drift@teststat["statistic", "tau2"] < qnorm(p=p.critical[critical])){
              # o coeficiente de y_t-1 ? diferente de zero, logo nao tem unit root
              dados.ret[i, "Unit.Root"] <- FALSE
              dados.ret[i, "Criterio"] <- 5
            } else {
              dados.ret[i, "Unit.Root"] <- TRUE
              dados.ret[i, "Criterio"] <- 6
            }
          }
        }
      } else {
        # aqui a trend ? significante, entao devemos estimar gamma com distrib normal
        if(ADF.trend@teststat["statistic", "tau3"] < qnorm(p=p.critical[critical])){
          dados.ret[i, "Unit.Root"] <- FALSE
          dados.ret[i, "Criterio"] <- 7
        }else{
          dados.ret[i, "Unit.Root"] <- TRUE  
          dados.ret[i, "Criterio"] <- 8
        }
      }
    }
    
    # Remove the models and variable x to avoid trash
    rm(list = c("x", "ADF.none", "ADF.drift", "ADF.trend"))
  }
  
  # Retorno dos dados
  return(dados.ret)
}


# User defined variables --------------------------------------------------


file.name <- "TCC BRUNO SANTOS-2008-2021- valores MENSAIS- atualiz. 24_05_2022.xlsx"


# Data load ---------------------------------------------------------------

SIH <- read_excel(file.name,
                  sheet = 1,
                  range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                  na = "-")

SIA <- read_excel(file.name,
                  sheet = 2,
                  range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                  na = "-")

SIH.Value <- read_excel(file.name,
                        sheet = 3,
                        range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                        na = "-")

SIA.Value <- read_excel(file.name,
                        sheet = 4,
                        range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                        na = "-")

IPCA <- read_excel(file.name,
                   sheet = "IPCA",
                   range = cell_limits(ul = c(2,1), lr = c(NA, NA)),
                   # col_types = c("date","numeric"),
                   na = "-")

IPCA$Data <- lubridate::ymd(IPCA$Data, truncated = 1)

# Data regularization -----------------------------------------------------

SIH$Data <- as.Date(SIH$Data)
SIA$Data <- as.Date(SIA$Data)
SIH.Value$Data <- as.Date(SIH.Value$Data)
SIA.Value$Data <- as.Date(SIA.Value$Data)


SIH %<>% select(c("Data", "Regiao_Norte", "Regiao_Nordeste", "Regiao_CentroOeste", "Regiao_Sudeste", "Regiao_Sul"))
SIA %<>% select(c("Data", "Regiao_Norte", "Regiao_Nordeste", "Regiao_CentroOeste", "Regiao_Sudeste", "Regiao_Sul"))
SIH.Value %<>% select(c("Data", "Regiao_Norte", "Regiao_Nordeste", "Regiao_CentroOeste", "Regiao_Sudeste", "Regiao_Sul"))
SIA.Value %<>% select(c("Data", "Regiao_Norte", "Regiao_Nordeste", "Regiao_CentroOeste", "Regiao_Sudeste", "Regiao_Sul"))


# SIH %<>% mutate(Data = as.Date(Data), DB = "SIH")
# SIA %<>% mutate(Data = as.Date(Data), DB = "SIA")

tbl_qtd <- SIA %>% bind_rows(SIH) %>% group_by(Data) %>% summarise_all(.funs = sum, na.rm=TRUE)

IPCA.last <- IPCA %>% filter(Data == "2021-01-01") %>% pull(IPCA)

# SIH.Value %>%
#   inner_join(IPCA, by = c("Data"="Data")) %>% 
#   mutate(fator = IPCA.last/IPCA,
#          Total_Brasil2 = Total_Brasil * fator )

tbl_value <- SIA.Value %>% 
  bind_rows(SIH.Value) %>%
  group_by(Data) %>% 
  summarise_all(.funs = sum, na.rm=TRUE) %>% 
  inner_join(IPCA, by = c("Data"="Data")) %>% 
  mutate(fator = IPCA.last/IPCA,
         Regiao_Norte = Regiao_Norte * fator, 
         Regiao_Nordeste = Regiao_Nordeste * fator, 
         Regiao_CentroOeste = Regiao_CentroOeste * fator,
         Regiao_Sudeste = Regiao_Sudeste * fator, 
         Regiao_Sul = Regiao_Sul * fator) %>% 
  select(-fator, - IPCA)


ggplot(tbl_qtd) + 
  geom_line(aes(x=Data, y=Regiao_Norte, colour = "Regiao_Norte")) +
  geom_line(aes(x=Data, y=Regiao_Nordeste, colour = "Regiao_Nordeste")) +
  geom_line(aes(x=Data, y=Regiao_CentroOeste, colour = "Regiao_CentroOeste")) +
  geom_line(aes(x=Data, y=Regiao_Sudeste, colour = "Regiao_Sudeste")) +
  geom_line(aes(x=Data, y=Regiao_Sul, colour = "Regiao_Sul")) +
  labs()

ggplot(tbl_value) + 
  geom_line(aes(x=Data, y=Regiao_Norte, colour = "Regiao_Norte")) +
  geom_line(aes(x=Data, y=Regiao_Nordeste, colour = "Regiao_Nordeste")) +
  geom_line(aes(x=Data, y=Regiao_CentroOeste, colour = "Regiao_CentroOeste")) +
  geom_line(aes(x=Data, y=Regiao_Sudeste, colour = "Regiao_Sudeste")) +
  geom_line(aes(x=Data, y=Regiao_Sul, colour = "Regiao_Sul")) +
  labs()


tbl_qtd %>%
  mutate(BR = Regiao_Norte + Regiao_Nordeste + Regiao_CentroOeste + Regiao_Sudeste + Regiao_Sul) %>% 
  Testa.RaizUnitaria(critical = "5pct")

tbl_value %>%
  mutate(BR = Regiao_Norte + Regiao_Nordeste + Regiao_CentroOeste + Regiao_Sudeste + Regiao_Sul) %>% 
  Testa.RaizUnitaria(critical = "5pct")

tbl_qtd %<>%
  mutate(BR = Regiao_Norte + Regiao_Nordeste + Regiao_CentroOeste + Regiao_Sudeste + Regiao_Sul,
         Pandemia1 = if_else(Data == as.Date("2020-03-01"), true = 1, false = 0), 
         Pandemia2 = if_else(Data >= as.Date("2020-03-01"), true = 1, false = 0))

tbl_value %<>%
  mutate(BR = Regiao_Norte + Regiao_Nordeste + Regiao_CentroOeste + Regiao_Sudeste + Regiao_Sul,
         Pandemia1 = if_else(Data == as.Date("2020-03-01"), true = 1, false = 0), 
         Pandemia2 = if_else(Data >= as.Date("2020-03-01"), true = 1, false = 0))


tbl_qtd %>%
  filter(Pandemia2 == 0) %>% 
  Testa.RaizUnitaria(critical = "1pct")

Y <- tbl_qtd %>% filter(Pandemia2 == 0) %>% pull(BR)
acf(Y, 24)
pacf(Y, 24)

for (i in 12:1) {
  mdl <- arima(Y,
               order = c(i, 0, 0))  
  
  assign(x = sprintf("mdl.%d", i), 
         value = mdl)
}


AIC(mdl.1, mdl.2, mdl.3, mdl.4, mdl.5, mdl.6, mdl.7, mdl.8, mdl.9, mdl.10, mdl.11, mdl.12)
BIC(mdl.1, mdl.2, mdl.3, mdl.4, mdl.5, mdl.6, mdl.7, mdl.8, mdl.9, mdl.10, mdl.11, mdl.12)


mdl.12$coef/diag(mdl.12$var.coef)^0.5

mdl <- arima(Y,
             order = c(12, 0, 0),
             fixed = c(NA, NA, NA, 0, 0, 0, NA, NA, 0, 0, NA, 0, NA))

abs(mdl$coef[colnames(mdl$var.coef)]/diag(mdl$var.coef)^0.5)

length(mdl.12$residuals)-(length(colnames(mdl.12$var.coef)) - length(colnames(mdl$var.coef)))
((sum(mdl$residuals^2) - sum(mdl.12$residuals^2))/((length(colnames(mdl.12$var.coef)) - length(colnames(mdl$var.coef)))))/(sum(mdl.12$residuals^2) /(length(mdl.12$residuals)-(length(colnames(mdl.12$var.coef))))) 

pred <- predict(mdl.12, n.ahead = 22)
pred <- predict(mdl, n.ahead = 22)

tbl_qtd$Prev <- NA
tbl_qtd$Prev[tbl_qtd$Pandemia2 == 1] <- pred$pred

tbl_qtd$PrevUpper <- NA
tbl_qtd$PrevUpper[tbl_qtd$Pandemia2 == 1] <- pred$pred + pred$se * 1

tbl_qtd$PrevLower <- NA
tbl_qtd$PrevLower[tbl_qtd$Pandemia2 == 1] <- pred$pred - pred$se * 1


ggplot(tbl_qtd) + 
  geom_line(aes(x=Data, y=BR)) + 
  geom_ribbon(aes(x=Data, ymax=PrevUpper, ymin = PrevLower),
              fill = "Blue",
              alpha = 0.3) +
  geom_line(aes(x=Data, y=Prev), colour = "Red", linetype = "solid") + 
  labs(title = "Procedimentos hospitalares do SUS",
       subtitle = "AIH aprovadas")
