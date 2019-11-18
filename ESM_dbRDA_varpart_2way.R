#The code was modified from the paper "Sokol, E.R., Herbold, C.W., Lee, C.K., Cary, S.C. and Barrett, J.E. (2013) Local and regional influences over soil microbial metacommunities in the Transantarctic Mountains. Ecosphere 4(11), art136."

require(vegan)
#-------------------------------------------------------------------------------
# -- given vars
dist.name<-"bray"
#dat.comm<-decostand(M,"hellinger")
dat.comm<-M
dat.pcnm<-S
dat.env<-E
#dat.nutrient=N
#-------------------------------------------------------------------------------
# -- get rid of potential confounding output vars
if(exists("dat.varpart")){rm(dat.varpart)}

# -- list of original vars in working environment to keep after analysis
inputvars<-ls()

#-------------------------------------------------------------------------------
# -- combined model, all vars
# -- note that scaling doesn't matter -- capscale probably scales predictor vars internally?
dat.env.pcnm<-as.data.frame(cbind(dat.env,dat.pcnm))

mod1.env.pcnm<-capscale(dat.comm~., distance=dist.name, dat.env.pcnm)

#-------------------------------------------------------------------------------
# -- fit environmental model, forward stepwise select variables
# clear values
varselect.env.list<-NULL
# --
mod0.env<-capscale(dat.comm~1, distance=dist.name, dat.env)
mod1.env<-capscale(dat.comm~., distance=dist.name, dat.env)

mod.step.env<-ordistep(mod0.env,scope=formula(mod1.env),direction="forward",perm.max=999)
varselect.env.list<-names(mod.step.env[[6]]$envcentre)

#-------------------------------------------------------------------------------
# -- fit spatial model, forward stepwise select variables
# clear values
varselect.pcnm.list<-NULL
# --
mod0.pcnm<-capscale(dat.comm~1, distance=dist.name, dat.pcnm)
mod1.pcnm<-capscale(dat.comm~., distance=dist.name, dat.pcnm)

mod.step.pcnm<-ordistep(mod0.pcnm,scope=formula(mod1.pcnm),direction="forward",perm.max=999)
varselect.pcnm.list<-names(mod.step.pcnm[[6]]$envcentre)

#-------------------------------------------------------------------------------
## -- combined model (env + pcnm) [abc]
# clear values
varselect.all.list<-NULL
# --
varselect.all.list<-c(varselect.env.list,varselect.pcnm.list)
form.env.AND.pcnm<-paste(ifelse(is.null(varselect.all.list),"~1","~1+"),
                         paste(varselect.all.list,collapse="+"))
mod.env.AND.pcnm<-update(mod1.env.pcnm,as.formula(form.env.AND.pcnm))

## -- env given pcnm model (env | pcnm) [a | c]
form.env.NOT.pcnm<-paste(
  paste(
    ifelse(is.null(varselect.env.list),"~1","~1+"),
    paste(varselect.env.list,collapse="+")),
  ifelse(is.null(varselect.pcnm.list),
         "",
         paste("+Condition(",paste(varselect.pcnm.list,collapse="+"),")")
         )
  )
mod.env.NOT.pcnm<-update(mod1.env.pcnm,as.formula(form.env.NOT.pcnm))

## -- pcnm given env model (pcnm | env) [c | a]
form.pcnm.NOT.env<-paste(
  paste(
    ifelse(is.null(varselect.pcnm.list),"~1","~1+"),
    paste(varselect.pcnm.list,collapse="+")),
  ifelse(is.null(varselect.env.list),
         "",
         paste("+Condition(",paste(varselect.env.list,collapse="+"),")")
         )
  )
mod.pcnm.NOT.env<-update(mod1.env.pcnm,as.formula(form.pcnm.NOT.env))

#-------------------------------------------------------------------------------
# -- fit model for env ~ space
form.envBYcommpcnm<-as.formula(
  paste(ifelse(is.null(varselect.pcnm.list),"dat.env~1","dat.env~1+"),
                         paste(varselect.pcnm.list,collapse="+"))
  )

mod.envBYcommpcnm<-capscale(form.envBYcommpcnm, distance="euclidean", dat.pcnm)
# mod0.envBYpcnm<-rda(form.envBYcommpcnm, dat.pcnm)
mod0.envBYpcnm<-capscale(scale(dat.env)~1, distance="euclidean", dat.pcnm)
mod1.envBYpcnm<-capscale(scale(dat.env)~., distance="euclidean", dat.pcnm)

mod.step.envBYpcnm<-ordistep(mod0.envBYpcnm,scope=formula(mod1.envBYpcnm),direction="forward",perm.max=999)
mod.step.envBYpcnm.anova<-anova(mod.step.envBYpcnm, step=1000,perm=1000)
varselect.envBYpcnm.list<-names(mod.step.envBYpcnm[[6]]$envcentre)



#-------------------------------------------------------------------------------
# -- Calculate varparts

dat.varpart<-data.frame(
  row.names=c("abc"),
  predictor.var=NA,adj.R2=NA,Df1=NA,Df2=NA,
  F=NA,N.Perm=NA,P.val=NA,
  vars.env=NA,vars.pcnm=NA,
  stringsAsFactors=FALSE)

# -- [abc] env AND pcnm
mod.temp<-mod.env.AND.pcnm
partition="abc"
predictor.var<-"env.AND.pcnm"
# --
if(!is.null(varselect.pcnm.list)||!is.null(varselect.env.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=RsquareAdj(mod.env.AND.pcnm)$adj.r.squared,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=ifelse(is.null(varselect.env.list),NA,paste(varselect.env.list,collapse=" + ")),
    vars.pcnm=ifelse(is.null(varselect.pcnm.list),NA,paste(varselect.pcnm.list,collapse=" + ")),
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=0,
    stringsAsFactors=FALSE)
}

# -- [ab] env
mod.temp<-mod.step.env
partition<-"ab"
predictor.var<-"env"
# --
if(!is.null(varselect.env.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=as.character(predictor.var),
    adj.R2=RsquareAdj(mod.temp)$adj.r.squared,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=ifelse(is.null(varselect.env.list),NA,paste(varselect.env.list,collapse=" + ")),
    vars.pcnm=NA,
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=0,
    stringsAsFactors=FALSE)
}

# -- [bc] pcnm
mod.temp<-mod.step.pcnm
partition<-"bc"
predictor.var<-"pcnm"
# --
if(!is.null(varselect.pcnm.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=as.character(predictor.var),
    adj.R2=RsquareAdj(mod.temp)$adj.r.squared,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=NA,
    vars.pcnm=ifelse(is.null(varselect.pcnm.list),NA,paste(varselect.pcnm.list,collapse=" + ")),
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=0,
    stringsAsFactors=FALSE)
}

# - calculate other partitions
abc<-dat.varpart["abc","adj.R2"]
ab<-dat.varpart["ab","adj.R2"]
bc<-dat.varpart["bc","adj.R2"]
b<-ab+bc-abc
a<-ab-b
c<-bc-b

# -- [a] env | pcnm
mod.temp<-mod.env.NOT.pcnm
partition<-"a"
predictor.var<-"env|pcnm"
# --
if(!is.null(varselect.pcnm.list)||!is.null(varselect.env.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=as.character(predictor.var),
    adj.R2=a,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=NA,
    vars.pcnm=NA,
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=a,
    stringsAsFactors=FALSE)
}

# -- [b] env intersect pcnm
dat.varpart["b","adj.R2"]<-b

# -- [c] pcnm | env
mod.temp<-mod.pcnm.NOT.env
partition<-"c"
predictor.var<-"pcnm|env"
# --
if(!is.null(varselect.pcnm.list)||!is.null(varselect.env.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=as.character(predictor.var),
    adj.R2=c,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=NA,
    vars.pcnm=NA,
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=c,
    stringsAsFactors=FALSE)
}

# -- [env.commpcnm] - env ~ pcnm selected for community data 
mod.temp<-mod.envBYcommpcnm
partition<-"env.commpcnm"
predictor.var<-"envBYcommpcnm"
# --
if(!is.null(varselect.pcnm.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=as.character(predictor.var),
    adj.R2=RsquareAdj(mod.temp)$adj.r.squared,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=NA,
    vars.pcnm=ifelse(is.null(varselect.pcnm.list),NA,paste(varselect.pcnm.list,collapse=" + ")),
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=0,
    stringsAsFactors=FALSE)
}

# -- [env] - env ~ pcnm selected for community and environmental data
mod.temp<-mod.step.envBYpcnm
partition<-"env"
predictor.var<-"envBYpcnm"
# --
if(!is.null(varselect.envBYpcnm.list)){
  anova.temp<-anova(mod.temp, step=1000,perm=1000)
  dat.varpart[partition,]<-data.frame(
    predictor.var=as.character(predictor.var),
    adj.R2=RsquareAdj(mod.temp)$adj.r.squared,
    Df1=anova.temp$Df[1],
    Df2=anova.temp$Df[2],
    F=anova.temp["Model","F"],
    N.Perm=substring(attr(anova.temp, "heading")[1], 93,95),
    P.val=anova.temp["Model","Pr(>F)"],
    vars.env=NA,
    vars.pcnm=ifelse(is.null(varselect.envBYpcnm.list),NA,paste(varselect.envBYpcnm.list,collapse=" + ")),
    stringsAsFactors=FALSE)
}else{
  dat.varpart[partition,c("predictor.var","adj.R2")]<-data.frame(
    predictor.var=predictor.var,
    adj.R2=0,
    stringsAsFactors=FALSE)
}

dat.varpart.mod.anova.env<-data.frame(mod.step.env$anova)
dat.varpart.mod.anova.pcnm<-data.frame(mod.step.pcnm$anova)

mod.list<-list(
  env=mod.step.env,
  pcnm=mod.step.pcnm)

# -- remove non-outputvars from global environment
rm(list=ls()[!ls()%in%c(inputvars,
                        "dat.varpart",
                        "dat.varpart.mod.anova.env",
                        "dat.varpart.mod.anova.pcnm",
                        "mod.list")])


