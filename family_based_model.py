def update_one_generation(pp,pw,ww,k,h,s):
    new_pp = (pp*pp*k*1*(1-s)
                +pw*pp*k*0.5*(1-s)
                +pp*pw*k*0.5*(1-h*s)
                +pw*pw*k*0.25*(1-h*s)
                +pp*(1-k)*1*(1-s)
                +pw*(1-k)*0.25*(1-h*s)
            )

    new_pw =( pw*pp*k*0.5*(1-s)
            +ww*pp*k*1*(1-s)
            +pp*pw*k*0.5*(1-h*s)
            +pw*pw*k*0.5*(1-h*s)
            +ww*pw*k*0.5*(1-h*s)
            +pp*ww*k*1*1
            +pw*ww*k*0.5*1
            +pw*(1-k)*0.5*(1-h*s)
            )

    new_ww = (ww*pw*k*0.5*(1-h*s)
            +ww*ww*k*1*1
            +ww*(1-k)*1*1
            )
    new_total = new_pp+new_pw+new_ww
    return new_pp/new_total,new_pw/new_total, new_ww/new_total


def calculate_allele_freq(pp,pw,ww):
    return pp+pw/2,ww+pw/2

def allele_to_genotype_freq(p):
    return p*p, 2*p*(1-p),(1-p)*(1-p)
