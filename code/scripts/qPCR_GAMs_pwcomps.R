# Pairwise comparisons


#This script computes the difference in smooth trends for each gene, and plots the comparisons


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


comp_VEGF <- difference_smooths(VEGF_gam, smooth = "s(WEEK)", newdata = newdat,
                                 unconditional = FALSE, frequentist = FALSE,
                                 n=100, partial_match = TRUE)

comp_HIF1 <- difference_smooths(HIF1_gam, smooth = "s(WEEK)", newdata = newdat,
                                 unconditional = FALSE, frequentist = FALSE,
                                 n=100, partial_match = TRUE)


comp_DEK <- difference_smooths(DEK_gam, smooth = "s(WEEK)", newdata = newdat,
                                unconditional = FALSE, frequentist = FALSE,
                                n=100, partial_match = TRUE)

comp_ALDOA <- difference_smooths(ALDOA_gam, smooth = "s(WEEK)", newdata = newdat,
                                  unconditional = FALSE, frequentist = FALSE,
                                  n=100, partial_match = TRUE)

comp_STAT3 <- difference_smooths(STAT3_gam, smooth = "s(WEEK)", newdata = newdat,
                                  unconditional = FALSE, frequentist = FALSE,
                                  n=100, partial_match = TRUE)


comp_PGK1 <- difference_smooths(PGK1_gam, smooth = "s(WEEK)", newdata = newdat,
                                 unconditional = FALSE, frequentist = FALSE,
                                 n=100, partial_match = TRUE)

comp_RPTOR <- difference_smooths(RPTOR_gam, smooth = "s(WEEK)", newdata = newdat,
                                  unconditional = FALSE, frequentist = FALSE,
                                  n=100, partial_match = TRUE)


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

#Plot all pairwise comparisons.


a1<-plot_pairwise_pointwise(comp_VEGF)+labs(y='VEGF')
b1<-plot_pairwise_pointwise(comp_HIF1)+labs(y='HIF-1')
c1<-plot_pairwise_pointwise(comp_DEK)+labs(y='DEK')
d1<-plot_pairwise_pointwise(comp_ALDOA)+labs(y='ALDOA')
e1<-plot_pairwise_pointwise(comp_STAT3)+labs(y='STAT3')
f1<-plot_pairwise_pointwise(comp_PGK1)+labs(y='PGK1')
g1<-plot_pairwise_pointwise(comp_RPTOR)+labs(y='RPTOR')

a1+b1+c1+d1+e1+f1+g1+plot_annotation(tag_levels = 'A')


# Significant differences exist for the following genes and pairs of groups:
#
#     - VEGF: CG-MTD, MET-MTD
# - HIF-1a: CG-MET, MET-MTD
# - ALDOA: CG-MET, MET-MTD
# - PGK1: CG-MET, CG-MTD, MET-MTD
# - RPTOR: MET-MTD, CG-MTD
#
# Total of 11 comparisons. In all cases there is only one interval where the differences exist.

#this function takes the model with the pairwise comparisons from 'difference_smooths', receives the input of which
#pairs of groups to use, and the condition
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
                axis.text=element_text(size=22))+
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
                axis.text=element_text(size=22))+
            thm1+
            facet_wrap(~Pair_group)

    }


}

p1<-si_sig_plot(comp_VEGF,'CG - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Vegf"))))
p2<-si_sig_plot(comp_VEGF,'MET - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Vegf"))))
p3<-si_sig_plot(comp_HIF1,'CG - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Hif-1a"))))
p4<-si_sig_plot(comp_HIF1,'MET - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Hif-1a"))))
p5<-si_sig_plot(comp_ALDOA,'CG - MET','upper_s<0')+ylab(expression(atop("Difference",paste("Aldoa"))))
p6<-si_sig_plot(comp_ALDOA,'MET - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Aldoa"))))
p7<-si_sig_plot(comp_RPTOR,'CG - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Raptor"))))
p8<-si_sig_plot(comp_RPTOR,'MET - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Raptor"))))
p9<-si_sig_plot(comp_PGK1,'CG - MET','upper_s<0')+ylab(expression(atop("Difference",paste("Pgk1"))))
p10<-si_sig_plot(comp_PGK1,'CG - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Pgk1"))))
p11<-si_sig_plot(comp_PGK1,'MET - MTD','lower_s>0')+ylab(expression(atop("Difference",paste("Pgk1"))))

#Now put all the plots together

# layout<-"
# AAABBBCCC
# AAABBBCCC
# AAABBBCCC
# DDDEEEFFF
# DDDEEEFFF
# DDDEEEFFF
# GGGHHHIII
# GGGHHHIII
# GGGHHHIII
# JJJ###KKK
# JJJ###KKK
# JJJ###KKK
# "

layout<-"
AAABBB###
AAABBBIII
AAABBBIII
CCCDDDIII
CCCDDDJJJ
CCCDDDJJJ
EEEFFFJJJ
EEEFFFKKK
EEEFFFKKK
GGGHHHKKK
GGGHHH###
GGGHHH###
"


#fix aesthetics for the plots

p1<-p1+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p2<-p2+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p3<-p3+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p4<-p4+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p5<-p5+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p6<-p6+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p7<-p7+labs(x="Week")+scale_x_continuous(breaks = br)
p8<-p8+labs(x="Week")+scale_x_continuous(breaks = br)
p9<-p9+labs(x="")+scale_x_continuous(labels = rep("", 3), breaks = br)
p11<-p11+labs(x="Week")+scale_x_continuous( breaks = br)

p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+plot_annotation(tag_levels = 'A')+plot_layout(design = layout)&
    theme(plot.tag = element_text(size=18,face='bold'),
          axis.title = element_text(size=14),
          axis.text = element_text(size=14))
