
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


tbl_qtd <- SIA %>%
  bind_rows(SIH) %>% 
  group_by(Data) %>%
  summarise_all(.funs = sum, na.rm=TRUE)

IPCA.last <- IPCA %>% 
  filter(Data == "2021-01-01") %>%
  pull(IPCA)

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
  # geom_line(aes(x=Data, y=Regiao_Norte, colour = "Regiao_Norte")) +
  # geom_line(aes(x=Data, y=Regiao_Nordeste, colour = "Regiao_Nordeste")) +
  # geom_line(aes(x=Data, y=Regiao_CentroOeste, colour = "Regiao_CentroOeste")) +
  # geom_line(aes(x=Data, y=Regiao_Sudeste, colour = "Regiao_Sudeste")) +
  # geom_line(aes(x=Data, y=Regiao_Sul, colour = "Regiao_Sul")) +
  geom_line(aes(x=Data, y=Total_Brasil, colour = "Total_Brasil")) +
  labs(title = "Quantidade")

ggplot(tbl_value) + 
  # geom_line(aes(x=Data, y=Regiao_Norte, colour = "Regiao_Norte")) +
  # geom_line(aes(x=Data, y=Regiao_Nordeste, colour = "Regiao_Nordeste")) +
  # geom_line(aes(x=Data, y=Regiao_CentroOeste, colour = "Regiao_CentroOeste")) +
  # geom_line(aes(x=Data, y=Regiao_Sudeste, colour = "Regiao_Sudeste")) +
  # geom_line(aes(x=Data, y=Regiao_Sul, colour = "Regiao_Sul")) +
  geom_line(aes(x=Data, y=Total_Brasil, colour = "Total_Brasil")) +
  labs(title = "Valor")

# Passa o log -------------------------------------------------------------

tbl_qtd <- tbl_qtd %>%
  mutate(ln.Br = log(Total_Brasil))

tbl_value <- tbl_value %>%
  mutate(ln.Br = log(Total_Brasil))

# Testa Raiz Unitaria -----------------------------------------------------

tbl <- tbl_qtd %>% 
  select(Data,
         Qtd = Total_Brasil,
         ln.Qtd = ln.Br) %>% 
  inner_join(tbl_value, by = c("Data")) %>% 
  select(Data,
         Qtd, ln.Qtd, 
         Value = Total_Brasil,
         ln.Value = ln.Br)

tbl %>% 
  Testa.RaizUnitaria(critical = "1pct")



# Adiciona dummies de pandemia --------------------------------------------

tbl %<>%
  mutate(Pandemia1 = if_else(Data == as.Date("2020-03-01"), true = 1, false = 0), 
         Pandemia2 = if_else(Data >= as.Date("2020-03-01"), true = 1, false = 0))


# Analise de efeito da pandemia no valor ----------------------------------
plot(tbl)
tbl <- tbl %>% mutate(trend = row_number())

ols <- lm(ln.Value ~ -1 + ln.Qtd*Pandemia2 + trend, data = tbl)

summary(ols)


ols <- lm(Value ~ -1 + Qtd * Pandemia2 , data = tbl)

summary(ols)

ols <- lm(ln.Value ~ -1 + ln.Qtd * Pandemia2 + Pandemia1, data = tbl)

summary(ols)



# Analise de serie temporal das Quantidades -------------------------------


Fteste <- function(model.u, model.r){
  
  # LR <- -2*(model.r$loglik - model.u$loglik)
  
  q <- dim(model.u$var.coef)[1] - dim(model.r$var.coef)[1]
  n <- length(model.u$residuals) - dim(model.u$var.coef)[1]
  
  f.stat <- ((sum(model.r$residuals^2) - sum(model.u$residuals^2))/q)/(sum(model.u$residuals^2)/n)
  return(pf(f.stat, q, n, lower.tail = FALSE))
}



tbl2 <- tbl %>% filter(Pandemia2 == 0)

mdl <- arima(tbl2$Qtd,
                order = c(0, 0, 0),
                # xreg = data.matrix(tbl[ , c("Pandemia1", "Pandemia2")]),
                # method = "ML",
                # fixed = c(NA, NA, NA,
                #           0, 0, 0,
                #           0, 0, NA,
                #           0, NA, NA,
                #           NA, NA, NA)
                transform.pars = TRUE)

acf(mdl$residuals, 36)
pacf(mdl$residuals, 36)
mdl
a <- abs(mdl$coef[colnames(mdl$var.coef)]/diag(mdl$var.coef)^0.5)
a
which.min(a)
# Fteste(model.u = mdl.12, model.r = mdl)
mdl




for (i in 15:1) {
  mdl <- arima(tbl2$Qtd,
               order = c(i, 0, 0),
               transform.pars = FALSE)  
  
  assign(x = sprintf("mdl.%d", i), 
         value = mdl)
}


tbl.aic <- AIC(mdl.1, mdl.2, mdl.3, mdl.4, mdl.5, mdl.6, mdl.7, mdl.8, mdl.9, mdl.10, mdl.11, mdl.12)
tbl.bic <- BIC(mdl.1, mdl.2, mdl.3, mdl.4, mdl.5, mdl.6, mdl.7, mdl.8, mdl.9, mdl.10, mdl.11, mdl.12)

tbl.aic[which.min(tbl.aic$AIC),]
tbl.bic[which.min(tbl.bic$BIC),]

abs(mdl.3$coef[colnames(mdl.3$var.coef)]/diag(mdl.3$var.coef)^0.5)

acf(mdl.11$residuals, 36)
pacf(mdl.11$residuals, 36)

abs(mdl.11$coef[colnames(mdl.11$var.coef)]/diag(mdl.11$var.coef)^0.5)


mdl.11r <- arima(tbl2$Qtd,
             order = c(11, 0, 0),
             fixed = c(NA, NA, NA,
                       0, 0, 0,
                       NA, NA, 0,
                       0, NA, NA),
             # method = "ML",
             transform.pars = FALSE)  
abs(mdl.11r$coef[colnames(mdl.11r$var.coef)]/diag(mdl.11r$var.coef)^0.5)

Fteste(model.u = mdl.11, model.r = mdl.11r)

mdl.11r

pred.1 <- predict(mdl.1, n.ahead = 22)
pred.3 <- predict(mdl.11r, n.ahead = 22)

mean(pred.3$pred)

tbl$Prev1 <- tbl$Qtd
tbl$Prev3 <- tbl$Qtd

tbl$Prev1[tbl$Pandemia2 == 1] <- pred.1$pred
tbl$Prev3[tbl$Pandemia2 == 1] <- pred.3$pred

tbl$Prev.1Upper <- NA
tbl$Prev.1Upper[tbl$Pandemia2 == 1] <- pred.1$pred + pred.1$se * 1

tbl$Prev.1Lower <- NA
tbl$Prev.1Lower[tbl$Pandemia2 == 1] <- pred.1$pred - pred.1$se * 1

tbl$Prev.3Upper <- NA
tbl$Prev.3Upper[tbl$Pandemia2 == 1] <- pred.3$pred + pred.3$se * 1.95

tbl$Prev.3Lower <- NA
tbl$Prev.3Lower[tbl$Pandemia2 == 1] <- pred.3$pred - pred.3$se * 1.95



ggplot(tbl) + 
  # geom_line(aes(x=Data, y=Qtd)) + 
  geom_ribbon(aes(x=Data, ymax=Prev.1Upper, ymin = Prev.1Lower),
              fill = "Blue",
              alpha = 0.3) +
  geom_line(aes(x=Data, y=Prev1), colour = "Red", linetype = "solid") + 
  labs(title = "Procedimentos hospitalares do SUS",
       subtitle = "AIH aprovadas")

ggplot(tbl) + 
  # geom_line(aes(x=Data, y=Qtd)) + 
  geom_ribbon(aes(x=Data, ymax=Prev.3Upper, ymin = Prev.3Lower),
              fill = "Blue",
              alpha = 0.3) +
  geom_line(aes(x=Data, y=Prev3), colour = "Red", linetype = "solid") + 
  labs(title = "Procedimentos hospitalares do SUS",
       subtitle = "AIH aprovadas")












# ANALISE DO VALOR  -------------------------------------------------------


for (i in 12:0) {
  mdl <- arima(tbl2$ln.Value,
               order = c(i, 0, 0),
               # method = "CSS",
               # SSinit = "Rossignol2011",
               transform.pars = FALSE, xreg = tbl2[, c("ln.Qtd")])
  
  assign(x = sprintf("mdl.%d", i), 
         value = mdl)
}


tbl.aic <- AIC(mdl.0, mdl.1, mdl.2, mdl.3, mdl.4, mdl.5, mdl.6, mdl.7, mdl.8, mdl.9, mdl.10, mdl.11, mdl.12)
tbl.bic <- BIC(mdl.0, mdl.1, mdl.2, mdl.3, mdl.4, mdl.5, mdl.6, mdl.7, mdl.8, mdl.9, mdl.10, mdl.11, mdl.12)

tbl.aic[which.min(tbl.aic$AIC),]
tbl.bic[which.min(tbl.bic$BIC),]

abs(mdl.1$coef[colnames(mdl.3$var.coef)]/diag(mdl.3$var.coef)^0.5)

acf(mdl.1$residuals, 36)
pacf(mdl.1$residuals, 36)

abs(mdl.11$coef[colnames(mdl.11$var.coef)]/diag(mdl.11$var.coef)^0.5)


mdl.11r <- arima(tbl2$ln.Value,
                 order = c(11, 0, 0),
                 fixed = c(NA, NA, NA,
                           0, 0, 0,
                           NA, NA, 0,
                           0, NA, NA),
                 # method = "ML",
                 transform.pars = FALSE)  
abs(mdl.11r$coef[colnames(mdl.11r$var.coef)]/diag(mdl.11r$var.coef)^0.5)

Fteste(model.u = mdl.11, model.r = mdl.11r)


pred.1 <- predict(mdl.1, n.ahead = 22)
pred.3 <- predict(mdl.11r, n.ahead = 22)

tbl$Prev1 <- tbl$ln.Value
tbl$Prev3 <- tbl$ln.Value

tbl$Prev1[tbl$Pandemia2 == 1] <- pred.1$pred
tbl$Prev3[tbl$Pandemia2 == 1] <- pred.3$pred

tbl$Prev.1Upper <- NA
tbl$Prev.1Upper[tbl$Pandemia2 == 1] <- pred.1$pred + pred.1$se * 1

tbl$Prev.1Lower <- NA
tbl$Prev.1Lower[tbl$Pandemia2 == 1] <- pred.1$pred - pred.1$se * 1

tbl$Prev.3Upper <- NA
tbl$Prev.3Upper[tbl$Pandemia2 == 1] <- pred.3$pred + pred.3$se * 1.95

tbl$Prev.3Lower <- NA
tbl$Prev.3Lower[tbl$Pandemia2 == 1] <- pred.3$pred - pred.3$se * 1.95



ggplot(tbl) + 
  # geom_line(aes(x=Data, y=Qtd)) + 
  geom_ribbon(aes(x=Data, ymax=Prev.1Upper, ymin = Prev.1Lower),
              fill = "Blue",
              alpha = 0.3) +
  geom_line(aes(x=Data, y=Prev1), colour = "Red", linetype = "solid") + 
  labs(title = "Procedimentos hospitalares do SUS",
       subtitle = "AIH aprovadas")

ggplot(tbl) + 
  # geom_line(aes(x=Data, y=Qtd)) + 
  geom_ribbon(aes(x=Data, ymax=Prev.3Upper, ymin = Prev.3Lower),
              fill = "Blue",
              alpha = 0.3) +
  geom_line(aes(x=Data, y=Prev3), colour = "Red", linetype = "solid") + 
  labs(title = "Procedimentos hospitalares do SUS",
       subtitle = "AIH aprovadas")




# Analise pos pandemia ----------------------------------------------------

mdl.11rf <- arima(tbl$ln.Qtd,
                 order = c(11, 0, 0),
                 fixed = c(NA, NA, NA,
                           0, 0, 0,
                           NA, NA, 0,
                           0, NA, NA, 0, NA),
                 xreg = tbl[, c("Pandemia1", "Pandemia2")],
                 # method = "ML",
                 transform.pars = FALSE)  

mdl.11rf

abs(mdl.11rf$coef[colnames(mdl.11rf$var.coef)]/diag(mdl.11rf$var.coef)^0.5)



mdl.11rf <- arima(tbl$ln.Value,
                  order = c(1, 0, 0),
                  # fixed = c(NA, NA, NA,
                  #           0, 0, 0,
                  #           NA, NA, 0,
                  #           0, NA, NA, 0, NA),
                  xreg = tbl[, c("ln.Qtd", "Pandemia1", "Pandemia2")],
                  # method = "ML",
                  transform.pars = FALSE)  

mdl.11rf
acf(mdl.11rf$residuals)
pacf(mdl.11rf$residuals)


