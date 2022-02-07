# qPCR Model diagnostics:

h<-appraise(VEGF_gam)+labs(title="VEGF")
i<-appraise(HIF1_gam)+labs(title="HIF-1a")
j<-appraise(DEK_gam)+labs(title="DEK")
k<-appraise(ALDOA_gam)+labs(title="ALDOA")
l<-appraise(STAT3_gam)+labs(title="STAT3")
m<-appraise(PGK1_gam)+labs(title="PGK1")
n<-appraise(RPTOR_gam)+labs(title="RPTOR")

(h|i)/(j|k)+plot_layout(byrow=T)

(l|m)/(n)+plot_layout(byrow=T)
