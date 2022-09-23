
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
                  sheet = "SIH- QUANT. FACO-FEC ",
                  range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                  na = "-")

SIA <- read_excel(file.name,
                  sheet = "SIA- QUANT APROVADA. FACO FEC",
                  range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                  na = "-")

SIH.Value <- read_excel(file.name,
                        sheet = "SIH- VALOR TOTAL FACO-FEC",
                        range = cell_limits(ul = c(17,1), lr = c(NA, NA)),
                        na = "-")

SIA.Value <- read_excel(file.name,
                        sheet = "SIA- VALOR APROVADO FACO-FEC",
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
  filter(Data == "2021-12-01") %>%
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

rm(list = c("file.name", "IPCA", "IPCA.last", "SIA", "SIA.Value", "SIH", 
            "SIH.Value", "tbl_qtd", "tbl_value"))


# Graficos ----------------------------------------------------------------


g1 <- ggplot(tbl) + 
  geom_line(aes(x=Data, y=Qtd, colour = "Atual")) +
  labs(title = "Procedimentos hospitalares e Produção Ambulatorial",
       subtitle = "SIH e SIA",
       x=NULL,
       y="Qtd. Aprovada",
       caption = "Fonte: Ministério da Saúde (SIA/SUS e SIH/SUS)", 
       colour = NULL) +
  theme_bw() +
  scale_color_manual(breaks = c("Atual", "Previsao", "Teo"),
                     values = c("#000000", "#FF0000", "#0000FF")) + 
  theme(legend.position = "bottom")

g1

ggsave(filename = "Atual SIA_SIH.png",
       plot = g1,
       units = "in",
       width = 8, height = 6,
       dpi = 100)


g2 <- ggplot(tbl) + 
  geom_line(aes(x=Data, y=Value, colour = "Atual")) +
  labs(title = "Valor financeiro aprovado referente a cirurgias de catarata",
       subtitle = "Nível de preço de Dezembro de 2021",
       x=NULL,
       y="Milhões [R$]",
       caption = "Fonte: Ministério da Saúde (SIA/SUS e SIH/SUS)", 
       colour = NULL) +
  theme_bw() +
  scale_color_manual(breaks = c("Atual", "Previsao", "Teo"),
                     values = c("#000000", "#FF0000", "#0000FF")) + 
  scale_x_date(breaks = seq(from = as.Date("2008-01-01"), to = as.Date("2022-01-01"), by="year"), 
               labels = 2008:2022) + 
  scale_y_continuous(breaks = seq(from = 1e7, to = 7e7, by=1e7), 
                     labels = sprintf("%d", seq(from = 1e7, to = 7e7, by=1e7)/1000000)) + 
  theme(legend.position = "none")

g2

ggsave(filename = "Valor financeiro cirurgia.png",
       plot = g2,
       units = "in",
       width = 8, height = 6,
       dpi = 100)



# Adiciona dummies de pandemia --------------------------------------------

tbl %<>%
  mutate(Pandemia1 = if_else(Data == as.Date("2020-03-01"), true = 1, false = 0), 
         Pandemia2 = if_else(Data >= as.Date("2020-03-01"), true = 1, false = 0))


summary(tbl)
sd(tbl$Qtd)
sd(tbl$Value)

# Analise de efeito da pandemia no valor ----------------------------------
tbl %>% select(Qtd, ln.Qtd, Value, ln.Value) %>% plot()

tbl <- tbl %>% mutate(Year.fc = factor(year(Data), levels = 2008:2021,
                                       labels = c("2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2009-2016",
                                                  "2017-2019",
                                                  "2017-2019",
                                                  "2017-2019",
                                                  "2020-2021",
                                                  "2020-2021")))

ggplot(tbl) + 
  geom_point(aes(x=Qtd, y=Value, colour = Year.fc)) + 
  # facet_wrap(~Year.fc, nrow = 3) +
  theme_bw() + 
  labs(colour=NULL) +
  theme(legend.position = "bottom")

tbl %>% 
  mutate(I1 = if_else(Year.fc == "2009-2016", true = 1, false = 0),
         I2 = if_else(Year.fc == "2017-2019", true = 1, false = 0),
         I3 = if_else(Year.fc == "2020-2021", true = 1, false = 0)) %>% 
  lm(Value ~ -1 + Qtd + I(Qtd*I2) + I(Qtd*I3), data = .) -> ols

summary(ols)




# Analise pos pandemia ----------------------------------------------------

mdl.11rf <- arima(tbl$Qtd,
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

length(mdl.11rf$residuals)

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

pred.11r <- predict(mdl.11r, n.ahead = 22)

mean(pred.11r$pred)

tbl$Prev1 <- tbl$Qtd
tbl$Prev3 <- tbl$Qtd

tbl$Prev1[tbl$Pandemia2 == 1] <- pred.1$pred
tbl$Prev3[tbl$Pandemia2 == 1] <- pred.11r$pred

tbl$Prev.1Upper <- NA
tbl$Prev.1Upper[tbl$Pandemia2 == 1] <- pred.1$pred + pred.1$se * 1

tbl$Prev.1Lower <- NA
tbl$Prev.1Lower[tbl$Pandemia2 == 1] <- pred.1$pred - pred.1$se * 1

tbl$Prev.3Upper <- NA
tbl$Prev.3Upper[tbl$Pandemia2 == 1] <- pred.11r$pred + pred.11r$se * 1.95

tbl$Prev.3Lower <- NA
tbl$Prev.3Lower[tbl$Pandemia2 == 1] <- pred.11r$pred - pred.11r$se * 1.95



ggplot(tbl) + 
  # geom_line(aes(x=Data, y=Qtd)) + 
  geom_ribbon(aes(x=Data, ymax=Prev.1Upper, ymin = Prev.1Lower),
              fill = "Blue",
              alpha = 0.3) +
  geom_line(aes(x=Data, y=Prev1), colour = "Red", linetype = "solid") + 
  labs(title = "Procedimentos hospitalares do SUS",
       subtitle = "AIH aprovadas")

range(tbl$Data)
g1 <- ggplot(tbl) + 
  geom_ribbon(aes(x=Data, ymax=Prev.3Upper, ymin = Prev.3Lower), fill = "Blue",
              alpha = 0.2) +
  geom_line(aes(x=Data, y=Prev3, colour = "Previsão"), linetype = "solid") + 
  geom_line(aes(x=Data, y=Qtd, colour = "Atual")) + 
  theme_bw() +
  labs(title = "Cirurgias de catarata",
       subtitle = "Previsão de SIH e SIA",
       x=NULL,
       y="Qtd. Aprovada",
       caption = "Fonte: Ministério da Saúde (SIA/SUS e SIH/SUS)", 
       colour = NULL) +
  scale_color_manual(breaks = c("Atual", "Previsão", "Teo"),
                     values = c("#000000", "#FF0000", "#0000FF")) + 
  scale_x_date(breaks = seq(from = as.Date("2008-01-01"), to = as.Date("2022-01-01"), by="year"), 
               labels = 2008:2022) + 
  theme(legend.position = "bottom")

g1

ggsave(filename = "Previsao SIA_SIH.png",
       plot = g1,
       units = "in",
       width = 8, height = 6,
       dpi = 100)


TotalPrevisto <- sum(tbl$Prev3[tbl$Pandemia2 == 1] )
TotalRealizado <- sum(tbl$Qtd[tbl$Pandemia2 == 1] )

TotalPrevisto - TotalRealizado

tbl %>% filter(Data >= "2020-03-01", Data <= "2020-12-01") %>% 
  pull(Prev3) %>% sum()

tbl %>% filter(Data >= "2020-03-01", Data <= "2020-12-01") %>% 
  pull(Qtd) %>% sum()


tbl %>% filter(Data >= "2020-03-01", Data <= "2020-12-01") %>% 
  pull(Value) %>% sum()

tbl %>% filter(Data >= "2020-03-01", Data <= "2020-12-01") %>% 
  pull(Prev3) %>% sum() * (614.68+187.51)

# --------------------

tbl %>% filter(Data >= "2020-03-01") %>% pull(Prev3) %>% sum()
tbl %>% filter(Data >= "2020-03-01") %>% pull(Prev3) %>% sum() * (614.68+187.51)



tbl %>% filter(Data >= "2020-03-01") %>% pull(Qtd) %>% sum()
tbl %>% filter(Data >= "2020-03-01") %>% pull(Value) %>% sum()







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




