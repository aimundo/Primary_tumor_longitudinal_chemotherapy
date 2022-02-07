
data1 <- read.csv(here("data","Extraction_merged_2018_2020_2021_cohorts_FINAL.csv"))

#change "DAY" to "WEEK"
data1<-rename(data1,WEEK=DAY)

data1$GROUP<-as.factor(data1$GROUP)
data1$HBO2<-(data1$HB*data1$STO2)
data1$HBO2_SD<-data1$HB_SD*data1$STO2
data1$HBO<-data1$HB*(1-data1$STO2)
data1$HBO_SD<-data1$HB_SD*(1-data1$STO2)

#change all MG labels to MET

data1<-data1 %>%
    mutate(GROUP=case_when(as.character(GROUP)=="MG"~"MET",
                           TRUE ~ as.character(GROUP)))

data1$GROUP<-as.factor(data1$GROUP)



#get mean values per group at wk1

mean_wk1 <- subset(data1, WEEK==1) %>%
    group_by(GROUP) %>%
    summarise(baseline_StO2=mean(STO2, na.rm=TRUE),
              baseline_HbO2=mean(HBO2, na.rm=TRUE),
              baseline_HbO=mean(HBO, na.rm=TRUE),
              baseline_tHb=mean(HB, na.rm=TRUE)
    )

#paste values back into the data1 dataframe
data1<-left_join(data1,mean_wk1,by="GROUP")

#calculate fold values

data1<-data1 %>%
    mutate(StO2_fold=STO2/baseline_StO2,
           tHb_fold=HB/baseline_tHb,
           HbO2_fold=HBO2/baseline_HbO2,
           HbO_fold=HBO/baseline_HbO)


# GAMs

#GAM with interaction of group and time. All models for fold-change values will follow the same basic construction:
#Gaussian case
# StO2_gam <- gam(StO2_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
#                 method='REML',
#                 data  = data1)
#
# HbO2_gam <- gam(HbO2_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
#                 method='REML',
#                 data  = data1)
#
# tHb_gam <- gam(tHb_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
#                method='REML',
#                data  = data1)
#
# HbO_gam <- gam(HbO_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
#                method='REML',
#                data  = data1)

    #model<-gam(Variable~Group+s(Day,by=Group=k=6,` \n `method= 'REML',`\n `data=data1)

#t-scaled case
StO2_gam <- gam(StO2_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
                family=scat(link='identity', min.df = 3),
                method='REML',
                data  = data1)

HbO2_gam <- gam(HbO2_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
                family=scat(link='identity', min.df = 3),
                method='REML',
                data  = data1)

tHb_gam <- gam(tHb_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
               family=scat(link='identity', min.df = 3),
               method='REML',
               data  = data1)

HbO_gam <- gam(HbO_fold ~ GROUP+s(WEEK, by = GROUP, k = 6),
               family=scat(link='identity', min.df = 3),
               method='REML',
               data  = data1)

# Plot GAMs


#breaks for plots
br<-seq(2,6,2)


GAM_plot2<-function(model,data){

    ci <- confint(model, parm = "s(WEEK)", partial_match = TRUE, type = "confidence")

    # simultaneous interval
    si <- confint(model, parm = "s(WEEK)", type = "simultaneous", partial_match = TRUE)

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
    si <- si %>%
    mutate(est = case_when(GROUP == "MET" ~ est + constMET,
                           GROUP == "MTD" ~ est + constMTD,
                                              TRUE ~ est+Intercept),
    lower = case_when(GROUP == "MET" ~ lower + constMET,
                      GROUP == "MTD" ~ lower + constMTD,
                                              TRUE ~ lower+Intercept),
    upper = case_when(GROUP == "MET" ~ upper + constMET,
                      GROUP == "MTD" ~ upper + constMTD,
                                              TRUE ~ upper+Intercept))


    plot<-ggplot(ci, aes(x = WEEK, y = est, group = smooth)) +
        # geom_violin(data=data,aes(x=WEEK,y=Y,group=WEEK,
        #                           color=GROUP),
        #             alpha=0.5,
        #             #color='gray',
        #             inherit.aes = FALSE,
        #             show.legend = FALSE)+
        geom_boxplot(data=data,aes(x=WEEK,y=Y, group=WEEK,
                                   color=GROUP),
                     alpha=0.5,
                     #color='gray',
                     inherit.aes = FALSE,
                     show.legend = FALSE)+
        geom_line(lwd = 1,show.legend=FALSE) +
        geom_ribbon(data = ci, mapping = aes(ymin = lower, ymax = upper, x = WEEK, group = smooth,fill = GROUP),
                    inherit.aes = FALSE, alpha = 0.8,
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
        #stat_summary(data=data,aes(x=WEEK,y=Y,group=GROUP),fun=(mean),geom="line")+
        #stat_summary(data=data,aes(x=WEEK,y=Y,group=GROUP),fun.data=mean_cl_normal,geom='errorbar',width=0.5)+
        geom_line(data=ci,aes(WEEK,upper,color=GROUP), size=0.8, alpha=0.7,show.legend = FALSE)+
        scale_x_continuous(breaks=br)+
        theme_classic()+
        theme(
            axis.text=element_text(size=22),
            strip.text=element_text(size=22))+
        facet_wrap(~GROUP)+
        thm1

}

a<-GAM_plot2(StO2_gam,data1)+labs(y=expression(StO[2]),x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
b<-GAM_plot2(HbO2_gam,data1)+labs(y=expression(atop(HbO[2])),x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
c<-GAM_plot2(tHb_gam,data1)+labs(y="tHb",x="Week")
d<-GAM_plot2(HbO_gam,data1)+labs(y="HbO",x="Week")

a+b+c+d+plot_annotation(tag_levels = "A")

#get estimates for each model and each group

#stO2
as.data.frame(a$data) %>% group_by(GROUP) %>%
    summarise(max=max(est,na.rm=TRUE))

#HbO2
as.data.frame(b$data) %>% group_by(GROUP) %>%
    summarise(max=max(est,na.rm=TRUE))

#tHb
as.data.frame(c$data) %>% group_by(GROUP) %>%
    summarise(max=max(est,na.rm=TRUE))

#HbO
as.data.frame(d$data) %>% group_by(GROUP) %>%
    summarise(max=max(est,na.rm=TRUE))


## Model diagnostics

#Plot model diagnostics (using {gratia}).

a1<-appraise(StO2_gam,ncol=4,nrow=1,point_col = 'steelblue',point_alpha = 0.5)+ggtitle("StO2")
b1<-appraise(HbO2_gam,ncol=4,nrow=1,point_col = 'steelblue',point_alpha = 0.5)+ggtitle("HbO2")
c1<-appraise(tHb_gam,ncol=4,nrow=1,point_col = 'steelblue',point_alpha = 0.5)+ggtitle("tHb")
d1<-appraise(HbO_gam,ncol=4,nrow=1,point_col = 'steelblue',point_alpha = 0.5)+ggtitle("HbO")

a1/b1/c1/d1+plot_annotation(tag_levels = 'A')+plot_layout(ncol=1,nrow=4,byrow = F)


## Pairwise comparisons


difference_pointwise <- function(f1, f2, smooth, by_var, smooth_var, data, Xp, V, coefs, nrep = 1000) {
    ## make sure f1 and f2 are characters
    f1 <-  as.character(f1)
    f2 <-  as.character(f2)
    cnames <- colnames(Xp)
    ## columns of Xp associated with pair of smooths
    c1 <- grepl(gratia:::mgcv_by_smooth_labels(smooth, by_var, f1), cnames, fixed = TRUE)
    c2 <- grepl(gratia:::mgcv_by_smooth_labels(smooth, by_var, f2), cnames, fixed = TRUE)
    ## rows of Xp associated with pair of smooths
    r1 <- data[[by_var]] == f1
    r2 <- data[[by_var]] == f2

    ## difference rows of Xp for pair of smooths
    X <- Xp[r1, ] - Xp[r2, ]

    ######IMPORTANT: uncommenting the following two lines
    #removes the group means from the comparison######

    ## zero the cols related to other splines
    # X[, ! (c1 | c2)] <- 0

    ## zero out the parametric cols
    #X[, !grepl('^s\\(', cnames)] <- 0

    ## compute difference
    sm_diff <- drop(X %*% coefs)
    se <- sqrt(rowSums((X %*% V) * X))
    nr <- NROW(X)

    ## Calculate posterior simulation for smooths
    coefs_sim <- t(rmvn(nrep, rep(0, nrow(V)), V))
    rownames(coefs_sim) <- rownames(V)
    simDev <- X %*% coefs_sim
    absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
    masd <- apply(absDev, 2, max)
    crit_s <- quantile(masd, prob = 0.95, type = 8)


    out <- list(smooth = rep(smooth, nr), by = rep(by_var, nr),
                level_1 = rep(f1, nr),
                level_2 = rep(f2, nr),
                diff = sm_diff, se = se,
                lower_s = sm_diff - crit_s * se,
                upper_s = sm_diff + crit_s*se)

    out <- new_tibble(out, nrow = NROW(X), class = "difference_smooth")
    ## Only need rows associated with one of the levels
    out <- bind_cols(out, data[r1, smooth_var])

    out
}

difference_smooths <- function(model,
                               smooth,
                               n = 100,
                               ci_level = 0.95,
                               newdata = NULL,
                               partial_match = TRUE,
                               unconditional = FALSE,
                               frequentist = FALSE,
                               nrep = 10000,
                               include_means = TRUE,
                               ...) {
    if (missing(smooth)) {
        stop("Must specify a smooth to difference via 'smooth'.")
    }

    # smooths in model
    S <- gratia::smooths(model) # vector of smooth labels - "s(x)"
    # select smooths
    select <-
        gratia:::check_user_select_smooths(smooths = S, select = smooth,
                                           partial_match = partial_match)#,
    # model_name = expr_label(substitute(object)))
    sm_ids <- which(select)
    smooths <- gratia::get_smooths_by_id(model, sm_ids)
    sm_data <- map(sm_ids, gratia:::smooth_data,
                   model = model, n = n, include_all = TRUE)
    sm_data <- bind_rows(sm_data)
    by_var <- by_variable(smooths[[1L]])
    smooth_var <- gratia:::smooth_variable(smooths[[1L]])
    pairs <- as_tibble(as.data.frame(t(combn(levels(sm_data[[by_var]]), 2)),
                                     stringsAsFactor = FALSE))
    names(pairs) <- paste0("f", 1:2)

    Xp <- predict(model, newdata = sm_data, type = "lpmatrix")
    V <- gratia:::get_vcov(model, unconditional = unconditional,
                           frequentist = frequentist)
    coefs <- coef(model)

    out <- pmap(pairs, difference_pointwise, smooth = smooth, by_var = by_var,
                smooth_var = smooth_var, data = sm_data, Xp = Xp, V = V,
                coefs = coefs, nrep = nrep)
    out <- bind_rows(out)
    crit <- qnorm((1 - ci_level) / 2, lower.tail = FALSE)

    out <- add_column(out,
                      lower = out$diff - (crit * out$se),
                      upper = out$diff + (crit * out$se),
                      .after = 6L)

    out <- out %>%
        mutate(Pair_group=as.factor(paste(level_1,"-",level_2)))
    out
}


#Now, get the comparisons.


comp_StO2 <- difference_smooths(StO2_gam, smooth = "s(WEEK)", newdata = newdat,
                                 unconditional = FALSE, frequentist = FALSE,
                                 n=100, partial_match = TRUE)


comp_tHb <- difference_smooths(tHb_gam, smooth = "s(WEEK)", newdata = newdat,
                                unconditional = FALSE, frequentist = FALSE,
                                n=100, partial_match = TRUE)

comp_HbO <- difference_smooths(HbO_gam, smooth = "s(WEEK)", newdata = newdat,
                                unconditional = FALSE, frequentist = FALSE,
                                n=100, partial_match = TRUE)


comp_HbO2 <- difference_smooths(HbO2_gam, smooth = "s(WEEK)", newdata = newdat,
                                 unconditional = FALSE, frequentist = FALSE,
                                 n=100, partial_match = TRUE)

#this function plots the comparisons

plot_pairwise_pointwise<-function(model){

    p11<-ggplot() +
        geom_line(data = model, aes(x = WEEK, y = diff),size=1, alpha=0.5) +
        geom_ribbon(data = model, aes(x = WEEK, ymin = lower_s, ymax = upper_s),
                    alpha = 0.5, inherit.aes = FALSE,
                    fill=rib_col,
                    show.legend=FALSE) +
        geom_hline(yintercept = 0, lty = 2, color = "red")+
        scale_x_continuous(breaks=br)+
        theme_classic()+
        theme(
            axis.text=element_text(size=22))+
        facet_wrap(~Pair_group)+
        thm1

    p11
}

a2<-plot_pairwise_pointwise(comp_StO2)+labs(y=expression(StO[2]))
b2<-plot_pairwise_pointwise(comp_HbO2)+labs(y=expression(HbO[2]))
c2<-plot_pairwise_pointwise(comp_tHb)+labs(y='tHb')
d2<-plot_pairwise_pointwise(comp_HbO)+labs(y='HbO')

a2+b2+c2+d2+plot_annotation(tag_levels = 'A')


#From the plots, the significant differences are only in: (Gaussian case)

#StO2: CG-MET
#HbO2: CG-MET
#tHb: CG-MET

#t-scaled case:

#StO2: CG-MET, CG-MTD
#HbO2: CG-MET
#tHb: CG-MET, MET-MTD

#this function gets the fitted object, pair and condition and plots it



si_sig_plot<-function(model,pair,cond){

    if (cond=='lower_s>0') {

        model2<-model %>%
            filter(Pair_group==pair) %>%
            filter(lower_s>0)

        init1=model2$WEEK[1]
        final1=model2$WEEK[nrow(model2)]

        p1<-model %>%
            filter(Pair_group==pair)%>%
            ggplot() +
            annotate("rect",
                     xmin =init1, xmax =final1,ymin=-Inf,ymax=Inf,
                     fill=rect_col,
                     alpha = 0.5)+
            geom_ribbon( aes(x = WEEK, ymin = lower_s, ymax = upper_s),
                         alpha = 0.8, inherit.aes = FALSE,
                         fill=rib_col,
                         show.legend=FALSE) +
            geom_line( aes(x = WEEK, y = diff, group=Pair_group),size=1, alpha=0.5) +
            geom_hline(yintercept = 0, lty = 2, color = "red")+
            scale_x_continuous(breaks=br)+
            theme_classic()+
            xlab("")+
            theme(
                axis.text=element_text(size=22),
                strip.text=element_text(size=22))+
            thm1+
            facet_wrap(~Pair_group)

    } else if (cond=='upper_s<0'){

        model2<-model %>%
            filter(Pair_group==pair) %>%
            filter(upper_s<0)

        init1=model2$WEEK[1]
        final1=model2$WEEK[nrow(model2)]

        plot<-model %>%
            filter(Pair_group==pair)%>%
            ggplot() +
            annotate("rect",
                     xmin =init1, xmax =final1,ymin=-Inf,ymax=Inf,
                     fill=rect_col,
                     alpha = 0.5)+
            geom_ribbon( aes(x = WEEK, ymin = lower_s, ymax = upper_s),
                         alpha = 0.8, inherit.aes = FALSE,
                         fill=rib_col,
                         show.legend=FALSE) +
            geom_line( aes(x = WEEK, y = diff, group=Pair_group),size=1, alpha=0.5) +
            geom_hline(yintercept = 0, lty = 2, color = "red")+
            scale_x_continuous(breaks=br)+
            theme_classic()+
            xlab("")+
            theme(
                axis.text=element_text(size=22),
                strip.text=element_text(size=22))+
            thm1+
            facet_wrap(~Pair_group)

    }


}

# layout<-"
# AAABBBEEE
# AAABBBEEE
# AAABBBFFF
# CCCDDDFFF
# CCCDDDGGG
# CCCDDDGGG
# "

# o1<-si_sig_plot(comp_StO2,'CG - MET','upper_s<0')+ylab(expression(atop("Difference",paste("StO"[2]))))+scale_x_continuous(labels = rep("", 3), breaks = c(2,4,6))
# o2<-si_sig_plot(comp_HbO2,'CG - MET','upper_s<0')+ylab(expression(atop("Difference",paste("HbO"[2]))))+scale_x_continuous(labels = rep("", 3), breaks = c(2,4,6))
# o3<-si_sig_plot(comp_tHb,'CG - MET','upper_s<0')+ylab("Difference\ntHb")

o1<-si_sig_plot(comp_StO2,'CG - MET','upper_s<0')+ylab(expression(atop("Difference",paste("StO"[2]))))+scale_x_continuous(labels = rep("", 3), breaks = br)+scale_y_continuous(breaks=seq(-0.15,0.05,0.1))
o2<-si_sig_plot(comp_StO2,'CG - MTD','upper_s<0')+ylab(expression(atop("Difference",paste("StO"[2]))))+scale_x_continuous(labels = rep("", 3), breaks = br)+scale_y_continuous(breaks=seq(-0.10,0.05,0.1))
o3<-si_sig_plot(comp_HbO2,'CG - MET','upper_s<0')+ylab(expression(atop("Difference",paste("HbO"[2]))))+scale_x_continuous(labels = rep("", 3), breaks = br)+scale_y_continuous(breaks=seq(-0.3,0.1,0.2))
o4<-si_sig_plot(comp_tHb,'CG - MET','upper_s<0')+ylab("Difference\ntHb")+scale_x_continuous(labels = rep("", 3), breaks = br)+scale_y_continuous(breaks=seq(-0.2,0.0,0.2))
o5<-si_sig_plot(comp_tHb,'MET - MTD','lower_s>0')+ylab("Difference\ntHb")+scale_x_continuous(breaks = br)+scale_y_continuous(breaks=seq(-0.1,0.3,0.1))+xlab("Week")

layout<-"
AAAABBBBEE
AAAABBBBEE
AAAABBBBFF
AAAABBBBFF
AAAABBBBGG
CCCCDDDDGG
CCCCDDDDHH
CCCCDDDDHH
CCCCDDDDII
CCCCDDDDII
"
# a+b+c+d+o1+o2+o3+o4+o5+plot_annotation(tag_levels = "A")+plot_layout(byrow = F,design=layout)&
#     theme(plot.tag = element_text(size=18,face='bold'),
#           axis.title = element_text(size=14),
#           axis.text = element_text(size=14))

