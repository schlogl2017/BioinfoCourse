library(vegan)
library(iNEXT)

data("BCI")

# qtas species por parcela
spc.arvores <- specnumber(BCI)

#media e desvio padrão
m <- mean(spc.arvores)
m
s <- sd(spc.arvores)
s

# grafico riqueza
# las = muda a rotação dos valores do eixo y de lateral
# para 'normal' ou horizontal
hist(spc.arvores, main='NUMERO DE ESPECIES', xlab = 'Numero de espécies', ylab = 'Frequência', las = 1, col = 'red')

# soma do num especies por lote (cada especie é uma coluna nos dados BCI)
num.ind.spc <- colSums(BCI)
num.ind.spc

max(num.ind.spc)
min(num.ind.spc)

# diagrama de ranking abundancia
# mostra a esturtura da distribuição das especies
num.ind.spc.ranking <- sort(num.ind.spc, decreasing = TRUE)

plot(num.ind.spc.ranking, type='l', las=1, xlab = 'Especies-sorted', ylab = 'Abundancia total')
# adicionando pontos
points(num.ind.spc.ranking, pch=1)

# quais especies com mais de 500 individuos
num.ind.spc.ranking[which(num.ind.spc.ranking > 500)]

# especies com menor numero
num.ind.spc.ranking[which(num.ind.spc.ranking == 1)]

# indices de diversidade (mais famosos Shannon)
shannon.index.BCI <- diversity(BCI, index = 'shannon')
mean(shannon.index.BCI)
sd(shannon.index.BCI)

# estimadores classicos de especies
specpool(BCI)

# curva de acumulo
acum.spc <- specaccum(BCI)
acum.spc
# graficos
plot(acum.spc)

# extrapolando a curva de acumulo
iNEXT(num.ind.spc)
plot(iNEXT(num.ind.spc))



































