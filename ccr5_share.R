#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CCR5 Wei & Nielsen replication
# Robert Maier
# 
# Analysis code for Maier, Akbari et al.
# "No statistical evidence for an effect of CCR5 - Delta 32 on lifespan in the UK Biobank cohort"
# DOI 10.1038/s41591-019-0710-1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(magrittr)
library(survival)
library(plink2R)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load genotype data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this assumes that there are three PLINK files with all the relevant SNPs

path = '~/Downloads/'

gt = read_plink(paste0(path, 'ukb_cal_chr3_rs62625034_2MB'))
gtimp = read_plink(paste0(path, 'plink_ccr5_3snps'))
gtexm = read_plink(paste0(path, 'ukb_fe_exm_ccr5_subs'))
snpsgt = c('rs62625034', 'rs113010081')
snpsimp = c('3:46414943_TACAGTCAGTATCAATTCTGGAAGAATTTCCAG_T', 'rs113010081')
snpsexm = c('3:46373452:D:32')

ukbphen31063 = read_table2(paste0(path, 'ukb31063.sample_qc.tsv.gz'))

wball = ukbphen31063 %>% filter(in.white.British.ancestry.subset == 1) %$% id 
gaukbbids = ukbphen31063 %>% filter(genotyping.array == 'UKBB') %$% id

rsdat = as.data.frame(gt$bed) %>%
  rownames_to_column() %>%
  set_colnames(c('s', gt$bim$V2)) %>%
  mutate(s=as.numeric(gsub(':.+','',s))) %>%
  select(s, snpsgt) %>%
  dplyr::rename(rs113010081_genotyped='rs113010081') %>%
  left_join(as.data.frame(gtimp$bed) %>%
              rownames_to_column() %>%
              set_colnames(c('s', gtimp$bim$V2)) %>%
              mutate(s=as.numeric(gsub(':.+','',s))) %>%
              select(s, snpsimp), by='s') %>%
  dplyr::rename(rs113010081_imputed='rs113010081') %>%
  left_join(as.data.frame(gtexm$bed) %>%
              rownames_to_column() %>%
              set_colnames(c('s', gtexm$bim$V2)) %>%
              mutate(s=as.numeric(gsub(':.+','',s))) %>%
              select(s, snpsexm), by='s') %>%
  dplyr::rename(rs333_imputed='3:46414943_TACAGTCAGTATCAATTCTGGAAGAATTTCCAG_T',
                rs62625034_genotyped='rs62625034')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load age of death data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

deathvars = c('age_at_recruitment'='21022',
              'date_of_attending_assessment_centre'='53',
              'year_of_birth'='34',
              'month_of_birth'='52',
              'age_at_death'='40007')

death_pheno = read_table2(paste0(path, 'ukb31063.raw_phenotypes_death.txt')) %>%
  set_colnames(c('s', names(deathvars)[match(names(.)[-1], paste0(deathvars, '-0.0'))])) %>%
  separate(date_of_attending_assessment_centre, c('recy', 'recm', 'recd'), '-', convert = T) %>%
  mutate(yob = year_of_birth+month_of_birth/12,
         recy2 = recy+recm/12,
         arec = recy2-yob,
         diff = arec-age_at_recruitment,
         deathyear = yob+age_at_death,
         agenow = max(deathyear, na.rm=T) - yob,
         agenow = ifelse(is.na(age_at_death), agenow, age_at_death))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival analysis - Ext Data Fig 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lvl2 = lvl = c('3:46373452:D:32',
               'rs62625034_genotyped',
               'rs113010081_genotyped',
               'rs113010081_imputed',
               'rs333_imputed')
lvl2[1] = 'rs333_sequenced'
names(lvl2) = lvl

death_pg = death_pheno %>%
  filter(s %in% wball) %>%
  transmute(s, death = (!is.na(age_at_death))+0, start=arec, end=pmax(start+0.01, agenow)) %>%
  inner_join(rsdat %>% gather(snp, gt, -s) %>% transmute(s, snp, gt), by='s')

coxpvals = death_pg %>% mutate(arr = 'all') %>%
  bind_rows(death_pg %>% filter(s %in% gaukbbids) %>% mutate(arr = 'new')) %>%
  group_by(arr, snp) %>%
  summarize(coxp = summary(coxph(Surv(time = start, time2 = end, event = death) ~ gt==2))$coefficients[1,5]/2) %>%
  ungroup %>% mutate(k = ' Cumulative survival rate', gt=NA, snp = ordered(snp, levels=lvl))

death_agegrp = death_pheno %>%
  select(s, arec, age_at_death, agenow) %>%
  slice(rep(1:nrow(.), 37)) %>%
  mutate(agegrp = rep(41:77, each = nrow(death_pheno))) %>%
  mutate(surv = (agegrp > arec & agegrp < (agenow-1))+0,
         died = agegrp==floor(age_at_death)+0,
         died = ifelse(is.na(died), 0, died)) %>%
  filter(s %in% wball) %>%
  inner_join(rsdat %>% gather(snp, gt, -s) %>% transmute(s, snp, gt), by='s') %>%
  mutate(agegrp = ifelse(agegrp == 77, 76, agegrp)) %>%
  group_by(agegrp, snp, gt) %>%
  summarize(surv=sum(surv), died=sum(died)) %>%
  group_by(snp, gt) %>%
  mutate(dr = died/surv, sr = 1-dr, cumsr = cumprod(sr)) %>%
  ungroup


#~~~~~~~~~~~~~~~~~~~~ Survival figure - Ext Data Fig 1 ~~~~~~~~~~~~~~~~~~~~

survival_figure = death_agegrp %>%
  filter(agegrp < 75 | snp != '3:46373452:D:32') %>%
  filter(snp != 'exomedel2') %>% mutate(` Cumulative survival rate`=cumsr) %>%
  mutate(` Survival rate` = 1-dr, `Death count` = died) %>%
  mutate(snp = ordered(snp, levels=lvl)) %>%
  gather(k, v, ` Cumulative survival rate`, `Death count`, ` Survival rate`) %>%
  mutate(gt=as.factor(gt)) %>%
  filter(k != 'Death count' | gt == 2)

survival_figure %>%
  ggplot(aes(agegrp, v, col=gt, group=gt)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(k ~ snp, scales='free', labeller = labeller(snp=lvl2), switch = 'y') +
  xlab('age') + ylab('') +
  geom_text(aes(label = round(coxp, 3)), x=60, y=0.85,
            data = coxpvals %>% filter(arr == 'all', snp != 'exomedel2')) +
  theme(panel.background = element_blank())


#~~~~~~~~~~~~~~~~~~~~ Sample counts - Supp Table 7 ~~~~~~~~~~~~~~~~~~~~

death_agegrp %>%
  filter(!is.na(gt), gt == 2, snp != 'exomedel2') %>%
  mutate(snp = ordered(snp, levels=lvl)) %>%
  select(agegrp, snp, died) %>%
  spread(snp, died) %>%
  View


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Concordance tables - Ext Data Fig 2, Supp Table 2, indirectly Supp Table 3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

concordance_data = rsdat %>%
  filter(s %in% wball) %>%
  dplyr::rename(exm=`3:46373452:D:32`) %>%
  filter(!is.na(exm)) %>%
  gather(snp, gt, -s, -exm)

#~~~~~~~~~~~~~~~~~~~~ Ext Data Fig 2 ~~~~~~~~~~~~~~~~~~~~

ext_data2 = concordance_data %>%
  filter(s %in% gaukbbids) %>%
  group_by(snp) %>% count(exm, gt) %>%
  group_by(snp, exm) %>% mutate(tot = sum(n)) %>%
  group_by(snp, exm, gt) %>% summarize(all = sum(n)/tot) %>%
  gather(cnt, n, all) %>% unite(cnt, cnt, exm) %>% spread(cnt, n, fill=0)

ext_data2 %>% View

#~~~~~~~~~~~~~~~~~~~~ Supp Table 2 ~~~~~~~~~~~~~~~~~~~~

supp_table2 = concordance_data %>%
  mutate(newarr = s %in% gaukbbids) %>%
  group_by(snp) %>% count(exm, gt, newarr) %>%
  group_by(snp, exm, gt) %>% summarize(new=n[newarr], all = sum(n)) %>%
  gather(cnt, n, all, new) %>% unite(cnt, cnt, exm) %>% spread(cnt, n, fill=0) %>%
  mutate(varrank = recode(snp,
                          rs62625034_genotyped=1,
                          rs113010081_genotyped=2,
                          rs113010081_imputed=3,
                          rs333_imputed=4)) %>%
  arrange(varrank, gt) %>% select(1:8)

supp_table2 %>% View


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HWE of linked variants - Ext Data Fig 3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hwe = function(n0, n1 = NULL, n2 = NULL, tab = FALSE) {
  # returns chi-square HWE p-value from genotype counts (or an observed/expected table)
  
  if(is.null(n1[1]) & is.null(n2[1])) {n1=n0[2]; n2=n0[3]; n0=n0[1]}
  n = n0+n1+n2
  chisq = n*((4*n0*n2-n1^2)/(2*n0+n1)/(2*n2+n1))^2
  if(!tab) return(pchisq(chisq, 1, lower.tail = FALSE))
  
  p = (n0+n1/2)/n
  e0 = n*p^2
  e1 = n*2*p*(1-p)
  e2 = n*(1-p)^2
  c2 = (n0-e0)^2/e0 + (n1-e1)^2/e1 + (n2-e2)^2/e2
  cbind(c(n0, n1, n2), c(e0, e1, e2))
}


n = 4e5
af = 0.112
r2 = 0.994
r2s = seq(0.7, 0.999, len=100)
hrs = c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
hrs = seq(0.7, 0.95, len=20)
reps = 10
reps = 100

out = tibble()
for(j in 1:length(hrs)) {
  hr = hrs[j]
  
  print(hr)
  
  for(i in 1:reps) {
    
    h1 = sample(0:1, n, replace=T, prob=c(1-af, af))
    h2 = sample(0:1, n, replace=T, prob=c(1-af, af))
    
    h3 = h1
    sel = runif(length(h3)) > r2
    h3[sel] = (h3[sel]+1)%%2
    h4 = h2
    sel = runif(length(h4)) > r2
    h4[sel] = (h4[sel]+1)%%2
    
    g1 = h1+h2
    g2 = h3+h4
    
    dropouts = which(g2==2)
    dropouts = dropouts[runif(length(dropouts)) > hr]
    r2out = cor(g1[-dropouts], g2[-dropouts], use='p')^2
    out %<>% bind_rows(tibble(hwe1 = hwe(c(table(g1[-dropouts]), 0)),
                              hwe2 = hwe(c(table(g2[-dropouts]), 0)), r2out, hr))
  }
}

#~~~~~~~~~~~~~~~~~~~~ HWE of linked variants - Ext Data Fig 3 ~~~~~~~~~~~~~~~~~~~~

out %>%
  rename(`SNP 2` = hwe1, `SNP 1` = hwe2) %>%
  gather(snp, p, `SNP 1`, `SNP 2`) %>%
  mutate(nlp = -log10(p)) %>%
  group_by(hr, snp) %>%
  summarize(p = mean(nlp), psd = sd(nlp), pl = quantile(nlp, 0.05), ph = quantile(nlp, 0.95)) %>%
  ggplot(aes(1-hr, p, col=snp, group=snp)) +
  geom_point() +
  geom_errorbar(aes(ymin=pl, ymax=ph), width=0) +
  xlab('Proportion of rare homozygous samples removed') +
  ylab('HWE -log10(p-value)') +
  theme(panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival power calculations - Ext Data Fig 4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reps = 100
#reps = 10
thresh = 0.05
inc = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
rrs = 1+inc
p_gt = rsdat %>% filter(s %in% wball, !is.na(rs62625034_genotyped)) %$% mean(rs62625034_genotyped == 2)
p_d = death_pheno %>% filter(s %in% wball) %$% mean(!is.na(age_at_death))

death_pg_rs62625034 = death_pg %>% filter(snp == 'rs62625034_genotyped') %>% select(-gt, -snp)

power = tibble()
set.seed(123)
for(rr in rrs) {
  print(rr)
  p_d_gt0 = p_d / (p_gt*rr + (1-p_gt))
  p_d_gt1 = p_d_gt0 * rr
  p_gt1_d = p_d_gt1 * p_gt / p_d
  p_gt1_a = (1-p_d_gt1) * p_gt / (1-p_d)
  
  for(i in 1:reps) {
    coxp = death_pg_rs62625034 %>%
      mutate(thisaf = ifelse(death == 0, p_gt1_a, p_gt1_d), gt = (runif(n()) < thisaf)+0) %>%
      summarize(p = summary(coxph(Surv(time=start, time2=end, event=death) ~ gt==1))$coefficients[1,5]/2) %$% p
    power %<>% bind_rows(tibble(rr, rep=i, coxp))
  }
}

powermn = power %>% group_by(rr) %>% summarize(mnsig = mean(coxp < thresh))

#~~~~~~~~~~~~~~~~~~~~ Power figure - Ext Data Fig 4 ~~~~~~~~~~~~~~~~~~~~

powermn %>%
  ggplot(aes(rr, mnsig)) +
  geom_point() +
  geom_line() +
  xlab('relative risk') +
  ylab('Power') +
  theme(panel.background = element_blank(),
        axis.line = element_line())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HWE - Supp Table 5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rsdat %>% mutate(samp = 'all') %>%
  bind_rows(rsdat %>% mutate(samp = 'exm') %>% filter(s %in% gtexm$fam$V1)) %>%
  filter(s %in% wball) %>%
  gather(snp, gt, -s, -samp) %>%
  count(snp, samp, gt) %>%
  filter(!is.na(gt)) %>%
  group_by(snp, samp) %>%
  summarize(p = hwe(n), n2obs = n[3], n2exp=((n[2]/2+n[3])/sum(n))^2*sum(n)) %>%
  ungroup %>%
  transmute(snp = ordered(snp, levels=lvl), samp,
            num = paste0(format.pval(p, eps=1e-60, digits=2), ' (', n2obs, ', ', round(n2exp), ')')) %>%
  spread(samp, num)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# phenotypic associations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# requires the following data frames with phenotypes:
# 'pcs', 'ukbphen31063' (covariates besides PCs), 'allphen', 'hormones', 'hormones', 'finn'
# also requires data frame 'phennames' with columns 'phenotype', 'variable_type', 'description'

phewas = rsdat %>%
  gather(k, gt, -s) %>%
  mutate(g01_2 = pmax(gt-1, 0)) %>%
  gather(k2, gt, -s, -k) %>% unite(snp, k, k2) %>% spread(snp, gt) %>%
  inner_join(pcs %>% select(s, isFemale, age, PC1:PC20), by='s') %>%
  inner_join(ukbphen31063 %>% transmute(s=id, genotyping.array), by='s') %>%
  inner_join(allphen, by=c('s'='userId')) %>%
  inner_join(hormones, by=c('s'='eid')) %>%
  inner_join(finn, by=c('s'='eid'))


snpcols = which(grepl('g01_2', names(phewas)))
covarcols = which(names(phewas) %in% c('genotyping.array', 'isFemale', 'age', paste0('PC', 1:20)))
phencols = (max(covarcols)+1):ncol(phewas)
bintraits = sapply(phewas %>% select(phencols), function(x) length(unique(na.omit(x))) == 2)
bintraitsnam = names(bintraits)[bintraits]
bincols = match(bintraitsnam, names(phewas))
contcols = setdiff(phencols, bincols)

covs = paste0(names(phewas)[covarcols], collapse=' + ')

for(i in snpcols) {
  out = tibble()
  snp = names(phewas)[i]
  covs2 = covs
  if(snp == 'rs113010081_genotyped_g01_2') covs2 = gsub(' \\+ genotyping.array', '', covs2)
  print(snp)
  for (j in bincols) {
    cat(paste('\r', length(bincols)-match(j, bincols)))
    phen = names(phewas)[j]
    if(j %% 10 == 0) print(length(phencols) - j)
    tryCatch({
      res = summary(glm(as.formula(paste0(phen, ' ~ `', snp, '` + ', covs2)),
                        data=phewas, family='binomial'))$coefficients[2, ]
      out %<>% bind_rows(data.frame(snp, phen, b=res[1], se=res[2], p=res[4], stringsAsFactors=F))
    }, error=function(e) NA)
  }
}

counts = tibble()
for(i in snpcols) {
  print(i)
  snp = names(phewas)[i]
  set.seed(123)
  d2 = phewas %>% filter(!!sym(snp) == 1) %>% select(bintraitsnam)
  d3 = phewas %>% filter(!!sym(snp) == 1) %>% select(contcols)
  counts %<>% bind_rows(data.frame(cnt=colSums(d2, na.rm=T)) %>% rownames_to_column(var='phen') %>% mutate(snp=snp))
  counts %<>% bind_rows(data.frame(cnt=colSums(!is.na(d3))) %>% rownames_to_column(var='phen') %>% mutate(snp=snp))
}

binmat = phewas %>% select(bintraitsnam) %>% as.matrix
cases = apply(binmat, 2, function(x) sum(x == max(x, na.rm=T), na.rm=T))
tot = apply(binmat, 2, function(x) sum(is.finite(x)))
prevdat = tibble(phen = names(cases), cases=cases, tot=tot) %>%
  transmute(phen, rate = cases/tot)

ext_data5 = out %>%
  left_join(prevdat, by='phen') %>%
  left_join(counts, by=c('snp', 'phen')) %>%
  mutate(phen = gsub('^X', '', phen)) %>%
  left_join(phennames %>% transmute(phen=gsub('_irnt', '', phenotype),
                                    variable_type, description), by='phen') %>%
  filter(!is.na(rate), cnt > 5) %>%
  group_by(snp) %>%
  mutate(or = exp(b), padj = p.adjust(p, method='bonferroni')) %>%
  arrange(snp, -padj) %>%
  ungroup

#~~~~~~~~~~~~~~ PheWAS - Ext Data Fig 5 ~~~~~~~~~~~~~~

ext_data5 %>%
  filter(b < 10) %>%
  mutate(snp = gsub('_g01_2', '', snp)) %>%
  ggplot(aes(rate, exp(b), col=-log10(p), shape = cnt > 10, phen=phen, description=description, cnt=cnt)) +
  geom_hline(yintercept=c(1, 1/4, 1/3, 1/2, 2, 3, 4), color=rep(c('black', rep('grey', 6)), 5)) +
  geom_point() +
  facet_wrap(~ snp, labeller=labeller(snp=lvl2), nrow=2) +
  viridis::scale_color_viridis() +
  scale_x_log10() +
  scale_y_log10() +
  xlab('Sample prevalence') +
  ylab('Odds ratio') +
  theme(panel.background = element_blank()) +
  scale_shape_manual(values=c(1,20))


# continuous phenotypes
covs2 = covs = paste0(names(phewas)[covarcols], collapse=' + ')
outlin = tibble()
for(i in snpcols) {
  snp = names(phewas)[i]
  print(snp)
  if(snp == 'rs113010081_genotyped_g01_2') covs2 = gsub(' \\+ genotyping.array', '', covs2)
  set.seed(123)
  for (j in contcols) {
    cat(paste('\r', length(contcols)-match(j, contcols)))
    phen = names(phewas)[j]
    tryCatch({
      res = summary(lm(as.formula(paste0('`', phen, '` ~ `', snp, '` + ', covs2)), data=phewas))$coefficients[2, ]
      outlin %<>% bind_rows(data.frame(snp, phen, b=res[1], se=res[2], p=res[4], stringsAsFactors=F))
    }, error=function(e) NA)
  }
}


#~~~~~~~~~~~~~~ QQ-plot - Ext Data Fig 6 ~~~~~~~~~~~~~~

get.lambda = function(p) {
  # calculate genomic inflaction factor lambda from p-values
  median(na.omit(qchisq(p, 1, lower.tail=F))) / .456
}

ggqq_df = function(dat, reduce=0, title = '', labs = TRUE, lambdagc = TRUE, confidence_region = TRUE) {
  # same as ggqq, but with (grouped) data.frame with column 'p' as input
  
  if(!is.null(groups(dat))) dat %<>%
    select(p, as.character(groups(.)[[1]])) %>%
    set_colnames(c('p', 'group')) %>%
    group_by(group)
  
  d2 = dat %>%
    filter(!is.na(p)) %>%
    transmute(p, pval = -log10(sort(p)),
              n=n(), a=1:n(), b=a/n(),
              upper=qbeta(0.025, a, rev(a)),
              lower = qbeta(0.975, a, rev(a)),
              x = -log10(b))
  
  if(reduce > 0) d2 %<>% mutate(sel = runif(n()) > reduce) %>% filter(x > 2 | sel)
  
  grp = ifelse('group' %in% names(d2), 'group', 'NULL')
  if(labs) {
    xl = substitute(paste(-log[10], "(expected P)"))
    yl = substitute(paste(-log[10], "(observed P)"))
  }
  
  pl = ggplot(d2, aes_string('x', 'pval', col=grp))
  if(confidence_region) pl = pl +
    geom_ribbon(aes_string(ymin='-log10(lower)',
                           ymax='-log10(upper)', fill=grp), alpha=.4, colour=NA)
  pl = pl +
    geom_point() +
    theme(plot.title = element_text(hjust = 0.5, size=8),
          panel.background = element_blank(),
          text = element_text(size=20),
          axis.line.x = element_line(),
          axis.line.y = element_line(color="black"),
          legend.title=element_blank()) +
    geom_abline(slope=1, intercept=0) +
    ggtitle(title) +
    xlab(xl) +
    ylab(yl)
  
  if(lambdagc) pl = pl + geom_text(data = summarize(d2, l=round(get.lambda(p), 2), y=max(pval)*0.9),
                                   aes(label=l, y=y), x=1.5, size=5)
  pl
}

out %>%
  mutate(type='binary') %>%
  bind_rows(outlin %>% mutate(type='continuous')) %>%
  left_join(counts %>% select(-perm), by=c('snp', 'phen')) %>%
  filter(is.na(cnt) | cnt > 5) %>%
  mutate(snp = gsub('_g01_2', '', snp),
         snp=ifelse(snp == '3:46373452:D:32', 'rs333_sequenced', snp)) %>%
  group_by(snp) %>%
  ggqq_df(lambdagc = F)


#~~~~~~~~~~~~~~ PheWAS results - Supp table 8 ~~~~~~~~~~~~~~
out %>%
  mutate(type='binary') %>%
  bind_rows(outlin %>% mutate(type='continuous')) %>%
  left_join(counts %>% select(-perm), by=c('snp', 'phen')) %>%
  filter(is.na(cnt) | cnt > 5) %>%
  mutate(snp = gsub('_g01_2', '', snp), snp=ifelse(snp == '3:46373452:D:32', 'rs333_sequenced', snp)) %>%
  mutate(phen = gsub('^X', '', phen)) %>%
  left_join(phennames %>% transmute(phen=gsub('_irnt', '', phenotype), variable_type, description), by='phen') %>%
  arrange(p) %>%
  filter(snp == 'rs113010081_genotyped') %>%
  View



#~~~~~~~~~~~~~~ Overall health rating - Supp table 9 ~~~~~~~~~~~~~~

tab = table(rs333$X2178, rs333$rs113010081_genotyped_g01_2)
expected = rowSums(tab)/sum(tab)*colSums(tab)[2]
cbind(tab, expected)


