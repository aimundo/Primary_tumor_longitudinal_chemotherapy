
#First load the data, remove outliers and missing observations, and calculate average Cq values
#for each sample. Also, 'GENERAL_GROUP' will become 'GROUP'.

data<-read.csv(here('data/All_samples_with_average_effs_for_all_plates_2021.csv'))
#remove rows where Cq=NaN
data<-data %>% filter(Cq!='NaN')
data<-data %>%filter(SAMPLE!=57) #to remove sample 57 that seems to be an outlier
data$GENE<-as.factor(data$GENE)
#remove 'GROUP' and rename 'GENERAL_GROUP' as GROUP
data<-select(data,-GROUP)

#rename column
data<-data %>%
    rename(GROUP=GENERAL_GROUP)


#keep only weeks 1, 4 and 6 for MTD samples

# data<-data %>%
#     filter(WEEK==c(1,4,6))
#change all MG labels to MET

data<-data %>%
    mutate(GROUP=case_when(as.character(GROUP)=="MG"~"MET",
                           TRUE ~ as.character(GROUP)))

data$GROUP<-as.factor(data$GROUP)

#Calculating average Cq per sample per gene
data1<-data %>%
    filter(GROUP %in% c('CG','MET','MTD')) %>%
    group_by(ID,GENE,WEEK) %>%
    mutate(mean_Cq=mean(Cq),Cq_SD=sd(Cq)) %>%
    filter(duplicated(mean_Cq)==FALSE)


#create dataframe and paste efficiencies
avg_effs<-tibble("GENE"=as.factor(c("VEGFA",
                                    "GAPDH",
                                    "HIF-1a",
                                    "RPTOR",
                                    "DEK",
                                    "ALDOA",
                                    "PGK1",
                                    "STAT3")),
                 "Eff"=c(1.97853947,
                         1.96252309,
                         2.0172954,
                         1.994233037,
                         1.98612065,
                         1.97435611,
                         1.91996324,
                         1.93909892)
)

#paste eff. values in the original dataframe
#using left_join: (x,y) ==> returns all rows from x and all columns from x and y.
data1<-data1 %>%
    left_join(avg_effs,data1,by="GENE")



#Sort the data1 dataframe alphabetically so it matches the order of the Av_Eff dataframe
data1<-data1[order(data1$GENE),]

#Calculate expression
data1 <-data1 %>%
    mutate(Expression=Eff^mean_Cq)


#Get the average Expression values for all genes for CONTROL WEEK 1,
#this will create a "calibrator" value, the mean expression per gene at week 1.
#Then, the expression per gene are extracted in new dataframes.


Calibrators<-data1 %>%
    group_by(GENE) %>%
    filter(WEEK==1,GROUP=='CG')%>%
    summarize(Baseline=mean(Expression))


#Values in new dataframes
#Data frames per gene
GAPDH<-data1 %>%
    filter(GENE %in% 'GAPDH')

VEGF<- data1 %>%
    filter(GENE %in% 'VEGFA')

HIF_1<- data1 %>%
    filter(GENE %in% 'HIF-1a')

RPTOR<- data1 %>%
    filter(GENE %in% 'RPTOR')

DEK<- data1 %>%
    filter(GENE %in% 'DEK')

PGK1<- data1 %>%
    filter(GENE %in% 'PGK1')

ALDOA<- data1 %>%
    filter(GENE %in% 'ALDOA')

STAT3<- data1 %>%
    filter(GENE %in% 'STAT3')

#Pasting the calibrator values in each dataframe

GAPDH<-GAPDH %>%
    mutate(Baseline=as.numeric(rep(Calibrators[3,2])))

VEGF<-VEGF%>%
    mutate(Baseline=as.numeric(rep(Calibrators[9,2])))

HIF_1<-HIF_1 %>%
    mutate(Baseline=as.numeric(rep(Calibrators[4,2])))

RPTOR<-RPTOR %>%
    mutate(Baseline=as.numeric(rep(Calibrators[6,2])))

DEK<-DEK %>%
    mutate(Baseline=as.numeric(rep(Calibrators[2,2])))

PGK1<-PGK1 %>%
    mutate(Baseline=as.numeric(rep(Calibrators[5,2])))

ALDOA<-ALDOA %>%
    mutate(Baseline=as.numeric(rep(Calibrators[1,2])))

STAT3<-STAT3 %>%
    mutate(Baseline=as.numeric(rep(Calibrators[7,2])))


#Pasting the GAPDH sample expression

VEGF$GAPDH_ex<-GAPDH$Expression
VEGF$GAPDH_Baseline<-GAPDH$Baseline

HIF_1$GAPDH_ex<-GAPDH$Expression
HIF_1$GAPDH_Baseline<-GAPDH$Baseline

RPTOR$GAPDH_ex<-GAPDH$Expression
RPTOR$GAPDH_Baseline<-GAPDH$Baseline

DEK$GAPDH_ex<-GAPDH$Expression
DEK$GAPDH_Baseline<-GAPDH$Baseline

PGK1$GAPDH_ex<-GAPDH$Expression
PGK1$GAPDH_Baseline<-GAPDH$Baseline

ALDOA$GAPDH_ex<-GAPDH$Expression
ALDOA$GAPDH_Baseline<-GAPDH$Baseline

STAT3$GAPDH_ex<-GAPDH$Expression
STAT3$GAPDH_Baseline<-GAPDH$Baseline

#Calculating Ratios

VEGF<-VEGF%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

HIF_1<-HIF_1%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

RPTOR<-RPTOR%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

DEK<-DEK%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

PGK1<-PGK1%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

ALDOA<-ALDOA%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))

STAT3<-STAT3%>%
    mutate(Ratio=(Baseline/Expression)/(GAPDH_Baseline/GAPDH_ex))



# With multiple samples, the mean moves beyond 1 and there is a need to
#normalize all ratio expression to the CG week 1 expression.
#In other words, all the ratios need to be normalized to CG week 1 ratio expression,
#which will make the mean of CG at week 1 effectively 1.
#First, all the gene dataframes are merged in a single dataframe and then extract
#the mean values of CG at week to normalize the values.


Norm<-bind_rows(VEGF,HIF_1,RPTOR,DEK,PGK1,ALDOA,STAT3)%>%
    group_by(GENE) %>%
    filter(WEEK==1,GROUP=='CG')%>%
    summarize(Norm_Ratio=mean(Ratio))

# Create the "Normalized ratio column by dividing each Ratio of expression by the CG week 1 values

VEGF<-VEGF%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[7])

HIF_1<-HIF_1%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[3])

RPTOR<-RPTOR%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[5])

DEK<-DEK%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[2])

PGK1<-PGK1%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[4])

ALDOA<-ALDOA%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[1])

STAT3<-STAT3%>%
    mutate(Norm_ratio=Ratio/Norm$Norm_Ratio[6])


#Fit the GAMs
k=6
VEGF_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
              method='REML',
              data=VEGF)

HIF1_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
              method='REML',
              data=HIF_1)

DEK_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
             method='REML',
             data=DEK)

ALDOA_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
               method='REML',
               data=ALDOA)

PGK1_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
              method='REML',
              data=PGK1)

RPTOR_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
               method='REML',
               data=RPTOR)

STAT3_gam<-gam(Norm_ratio~GROUP+s(WEEK,by=GROUP,k=k),
               method='REML',
               data=STAT3)


#Now, use a function to fit a GAM and plot the model and the data


#breaks for plots
br<-seq(2,6,2)

GAM_plot2<-function(model,data){

    ci <- confint(model, parm = "s(WEEK)", partial_match = TRUE, type = "confidence")

    ## simultaneous interval
    #si <- confint(model, parm = "s(WEEK)", type = "simultaneous", partial_match = TRUE)

    # mean shift for Treatment group
    Intercept <- coef(model)[1]
    constMET <- coef(model)[2]+Intercept
    constMTD <- coef(model)[3]+Intercept

    #extract the response from the model
    test<-model$formula[[2]]

    #convert to character
    test<-as.character(test)
    test<-as.list(test)
    #we are doing fold changes, so split the string where the diagonal is. We will end up
    #with a list of lenght three, only the last two contain the terms of interest


    #select columns that match the pattern
    data<-data %>%
        select(WEEK,GROUP,"Y"=test[[1]])


    #pointwise confidence interval
    ci <- ci %>%
        mutate(est = case_when(GROUP == "MET" ~ est + constMET,
                               GROUP == "MTD" ~ est + constMTD,
                               TRUE ~ est+Intercept),
               lower = case_when(GROUP == "MET" ~ lower + constMET,
                                 GROUP == "MTD" ~ lower + constMTD,
                                 TRUE ~ lower+Intercept),
               upper = case_when(GROUP == "MET" ~ upper + constMET,
                                 GROUP == "MTD" ~ upper + constMTD,
                                 TRUE ~ upper+Intercept))

    #simultaneous confidence interval
    # si <- si %>%
    # mutate(est = case_when(GROUP == "MET" ~ est + constMET,
    #                        GROUP == "MTD" ~ est + constMTD,
    #                                           TRUE ~ est+Intercept),
    # lower = case_when(GROUP == "MET" ~ lower + constMET,
    #                   GROUP == "MTD" ~ lower + constMTD,
    #                                           TRUE ~ lower+Intercept),
    # upper = case_when(GROUP == "MET" ~ upper + constMET,
    #                   GROUP == "MTD" ~ upper + constMTD,
    #                                           TRUE ~ upper+Intercept))
    #

    GAM_plot<-ggplot(ci, aes(x = WEEK, y = est, group = smooth)) +
        geom_boxplot(data=data, aes(x = WEEK,
                                    y = Y,
                                    color = GROUP,
                                    group=WEEK),
                     alpha=0.5,
                     inherit.aes = FALSE,
                     show.legend = FALSE)+
        geom_line(lwd = 1,show.legend=FALSE) +
        geom_ribbon(data = ci, mapping = aes(ymin = lower, ymax = upper, x = WEEK, group = smooth,fill = GROUP),
                    inherit.aes = FALSE, alpha = 0.7,
                    show.legend=FALSE) +
        # geom_ribbon(data = si,
        # mapping = aes(ymin = lower, ymax = upper, x = WEEK, group = smooth,fill =GROUP),
        # inherit.aes = FALSE, alpha = 0.3,
        # show.legend=FALSE)+
        # geom_point(data=data, aes(x = WEEK,
        #                           y = Y,
        #                           color = GROUP),
        #            size=1.5,
        #            alpha=0.6,
        #            inherit.aes = FALSE,
        #            show.legend = FALSE)+
        geom_line(data=ci,aes(WEEK,upper,color=GROUP), size=0.8, alpha=0.7,show.legend = FALSE)+
        scale_x_continuous(breaks=br)+
        theme_classic()+
        labs(y="",x="")+
        theme(
            axis.text=element_text(size=22),
            strip.text=element_text(size=22))+
        facet_wrap(~GROUP)+
        thm1

}




#Plotting all models:

layout<-"
AAABBBCCC
AAABBBCCC
AAABBBCCC
DDDEEEFFF
DDDEEEFFF
DDDEEEFFF
###GGG###
###GGG###
###GGG###
"
a<-GAM_plot2(VEGF_gam,VEGF)+labs(title="Vegf",y='Relative\nExpression')+scale_x_continuous(labels = rep("", 3), breaks = c(2,4,6))
b<-GAM_plot2(HIF1_gam,HIF_1)+labs(title="Hif-1a")+scale_x_continuous(labels = rep("", 3), breaks = c(2,4,6))
c<-GAM_plot2(DEK_gam,DEK)+labs(title="Dek")+scale_x_continuous(labels = rep("", 3), breaks = c(2,4,6))
d<-GAM_plot2(ALDOA_gam,ALDOA)+labs(title="Aldoa",y='Relative\nExpression',x='Week')
e<-GAM_plot2(STAT3_gam,STAT3)+labs(title="Stat3")+scale_x_continuous(labels = rep("", 3), breaks = c(2,4,6))
f<-GAM_plot2(PGK1_gam,PGK1)+labs(title="Pgk1",x='Week')
g<-GAM_plot2(RPTOR_gam,RPTOR)+labs(title="Raptor",y='Relative\nExpression',x='Week')

# a+b+c+d+e+f+g+plot_annotation(tag_levels = "A")+plot_layout(design = layout)&
#     theme(plot.tag = element_text(size=18,face='bold'))&
#     theme(axis.title = element_text(size=18))

#get values at week 1 for all groups

#VEGF
(VEGF) %>%
    filter(WEEK==1)%>%
    group_by(GROUP) %>%
    summarise(mn=mean(Norm_ratio))

#HIF-1a
(HIF_1) %>%
    filter(WEEK==1)%>%
    group_by(GROUP) %>%
    summarise(mn=mean(Norm_ratio))

