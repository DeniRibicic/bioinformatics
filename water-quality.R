library (tidyverse)
library (ggplot2)
library (viridis)
library (lubridate)
library (reshape2)
########################
setwd("/home/denir/Projects/MonMic/water_quality/")
plot.new()
#Plot
df <- read.table('dfm_belsvik_ST21.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df$value <- as.numeric(df$value)
#+ fig.height=6, fig.width=10
p <- ggplot(data=df, aes(x=as.Date(Dato), y=value,color=variable))+
  geom_point()+
  geom_smooth(method = 'lm')+
  scale_color_viridis(discrete = TRUE)+
  facet_wrap(~variable, scales="free_y")+
  xlab("")+
  ggtitle('B_21')+
  theme(legend.position="none") +
  scale_fill_viridis()
p +theme(axis.text.x=element_text(angle=45,hjust=1),
         strip.text = element_text(size=9))

df <- read.table('dfm_belsvik_ST25.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df$value <- as.numeric(df$value)

#+ fig.height=6, fig.width=10
p <- ggplot(data=df, aes(x=as.Date(Dato), y=value,color=variable))+
  geom_point()+
  geom_smooth(method = 'lm')+
  scale_color_viridis(discrete = TRUE)+
  facet_wrap(~variable, scales="free_y")+
  xlab("")+
  ggtitle('B_25')+
  theme(legend.position="none") +
  scale_fill_viridis()
p +theme(axis.text.x=element_text(angle=45,hjust=1),
         strip.text = element_text(size=9))


# Bremnnes
df <- read.table('dfm_bremnes.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df$value <- as.numeric(df$value)
df$value <- as.numeric(df$value)

#+ fig.height=6, fig.width=10
p <- ggplot(data=df, aes(x=as.Date(Dato), y=value,color=variable))+
  geom_point()+
  geom_smooth(method = 'lm')+
  scale_color_viridis(discrete = TRUE)+
  facet_wrap(~variable, scales="free_y")+
  xlab("")+
  ggtitle('C')+
  theme(legend.position="none") +
  scale_fill_viridis()
p +theme(axis.text.x=element_text(angle=45,hjust=1),
         strip.text = element_text(size=9))



#Laksefjord
df <- read.table('dfm_laksefjord.csv',header=T,sep=',',dec='.',fileEncoding="UTF-8",na.strings = '')
df = na.omit(df)
df = subset(df,df$variable!='Dag')
df$value <- as.numeric(df$value)

#+ fig.height=6, fig.width=10
p <- ggplot(data=df, aes(x=as.Date(Dato), y=value,color=variable))+
  geom_point()+
  geom_smooth(method = 'lm')+
  scale_color_viridis(discrete = TRUE)+
  facet_wrap(~variable, scales="free_y")+
  xlab("")+
  ggtitle('A')+
  theme(legend.position="none") +
  scale_fill_viridis()
p +theme(axis.text.x=element_text(angle=45,hjust=1),
         strip.text = element_text(size=10))

