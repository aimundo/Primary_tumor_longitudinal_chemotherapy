
#read data
data<-read.csv(here("data","Imaging_data.csv"))
data$Group<-as.factor(data$Group)

jitter <- position_jitter(width = 0.1, height = 0.1,seed=2021)

a1<-ggplot(data=data,aes(x=Group,y=MVD,fill=Group,color=Group))+
    geom_violin(alpha=0.5,show.legend=FALSE,width=1)+
    geom_boxplot(width=.1,color='black',outlier.shape = NA,size=1,
                 alpha=0.5,
                 show.legend = FALSE)+
    geom_point(show.legend=FALSE,shape=21,colour='black',size=2,position=jitter,
               alpha=0.5)+
    theme_classic()+
    geom_signif(
        comparisons=list(c("MET","MTD")),
        y_position = 43,
        textsize = 3,
        color='black',
        test='wilcox.test',
        map_signif_level = TRUE
    )+
    geom_signif(
        comparisons=list(c("CG","MTD")),
        y_position = 48,
        textsize = 3,
        test='wilcox.test',
        color='black',
        map_signif_level = TRUE
    )+
    scale_x_discrete(expand=c(0,0))+
    labs(y="MVD(%)")+
    ylim(0,55)+
    thm1

#Linear mixed model for MVD (a one-way ANOVA?) with Bonferroni correction

#mod<-lm(MVD~Group,data=data)
#post_hoc<-pairs(emmeans(mod,"Group",adjust="bonf"))

#Kruskal Wallis test with Bonferroni correction as well

mod1<-kruskal.test(MVD~Group,data=data)
pairwise.wilcox.test(data$MVD,data$Group,p.adjust.method = "bonferroni")


###get median values of NEstin

data %>% filter(Group=='MTD') %>% summarise(median=median(MVD))
data %>% filter(Group=='CG') %>% summarise(median=median(MVD))
data %>% filter(Group=='MET') %>% summarise(median=median(MVD))

## add pictures to create composite image
CG_Nestin<-rasterGrob(readPNG(here("figures","CG_Nestin-scaled.png"),native=TRUE))
CG_MET<-rasterGrob(readPNG(here("figures","MET_Nestin-scaled.png"),native=TRUE))
CG_MTD<-rasterGrob(readPNG(here("figures","MTD_Nestin-scaled.png"),native=TRUE))

image_a<-wrap_elements(
    panel=CG_Nestin
)+
    ggtitle('CG')+
    theme(plot.title = element_text(hjust = 0.5))+
    wrap_elements(
        panel=CG_MET
    )+
    ggtitle('MET')+
    theme(plot.title = element_text(hjust = 0.5))+
    wrap_elements(
        panel=CG_MTD
    )+
    ggtitle('MTD')+
    theme(plot.title = element_text(hjust = 0.5))+
    a1+
    plot_annotation(tag_levels = 'A')
