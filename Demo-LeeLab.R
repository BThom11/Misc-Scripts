library(ggplot2)
library(dplyr)
library(tidyr)
Info=read.csv("C:/Users/rjtho/Downloads/September-29-November-6-2016-Information-engaged-wary/September 29 - November 6, 2016 - Information engaged wary/InfoEngage.csv")

#Age and sex
ggplot(Info,aes(x=age,fill=factor(Info$sex,labels=c("Male","Female"))))+
  geom_bar()+
  theme(legend.title=element_blank())+
  xlab('Age (98=Dont know, 99=refuse)')

#Q10, How do you make major life decisions
Q10=Info[,c(1,68:74)]
Q10= Q10 %>% gather(question,ans,colnames(Q10[,2:8]))
xlabs=c("Stick with it","Follow Gut",'If Doubts, I retrace steps','Communicate to Others','Consider variety of situations','Weight Important Factors','Dont need background research')

ggplot(Q10,aes(x=factor(question,labels=xlabs),fill=factor(ans,labels=c("yes",'no','dont know','refuse'))))+
  geom_bar()+
  theme(legend.title=element_blank())+
  ggtitle('How do you make major life decisions?')

#Q9 How frequently do you feel you cant gather enough info?
ggplot(Info,aes(x=factor(q9,labels=c('Freq','Sometimes','Not Too Often','Never','Dont Know','Refuse'))))+
  geom_bar()+
  xlab(element_blank())+
  ggtitle("How Frequently do you feel you can NOT gather enough info for important decisions")

#Q6 How much do you trust following info sources
Q6=Info[,c(which(names(Info)=='q6a'):which(names(Info)=='q6h'))]
Q6= Q6 %>% gather(q,a,colnames(Q6))
xlabq6=c('National News','Social Media','Family/Friends','From Libraries','Local News','Government','Health Care Providers','Financial Inst.')

ggplot(Q6,aes(x=factor(q,labels=xlabq6),fill=factor(a,labels=c('Alot','Some','Not Much','Not at all','DK','REF'))))+
  geom_bar()+
  theme(legend.title=element_blank())+
  xlab(element_blank())+
  ggtitle("How much do you trust these Info sources")
