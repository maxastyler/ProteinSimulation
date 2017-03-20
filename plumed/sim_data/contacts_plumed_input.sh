cgcontacts=$1
cat $cgcontacts | grep -v "#" | awk '
 BEGIN{
  for(i=1;i<=9;i++)  ntail[i]=1
  for(i=10;i<=26;i++) alpha[i]=1
  for(i=27;i<=30;i++) coil[i]=1
  for(i=31;i<=62;i++) beta[i]=1
  for(i=63;i<=81;i++) coil[i]=1
  for(i=82;i<=106;i++) beta[i]=1
  for(i=107;i<=107;i++) ctail[i]=1
 }
 {
  n++
  t1=$1
  r1=$2
  i1=$3
  t2=$4
  r2=$5
  i2=$6
  d=$7
  w=$8
  all[n]=w
  # print "D"n": DISTANCE ATOMS="i1","i2
  # plumed works in nm!
  print "C"n": COORDINATION GROUPA="i1","i2" R_0="0.12*d" NN=8 MM=10"
  # ntail_glob
  if(((r1 in ntail)&&((r2 in ntail)||(r2 in alpha)||(r2 in coil)||(r2 in beta)))||((r2 in ntail)&&((r1 in ntail)||(r1 in alpha)||(r1 in coil)||(r1 in beta)))) ntail_glob[n]=w
  # ctail_glob
  if(((r1 in ctail)&&((r2 in ctail)||(r2 in alpha)||(r2 in coil)||(r2 in beta)))||((r2 in ctail)&&((r1 in ctail)||(r1 in alpha)||(r1 in coil)||(r1 in beta)))) ctail_glob[n]=w
  # ntail_ctail
  if(((r1 in ntail)&&(r2 in ctail))||((r2 in ntail)&&(r1 in ctail))) ntail_ctail[n]=w
  # alpha_alpha
  if(((r1 in alpha)&&(r2 in alpha))) alpha_alpha[n]=w
  # beta_beta
  if(((r1 in beta)&&(r2 in beta))) beta_beta[n]=w
  # alpha_beta
  if(((r1 in alpha)&&(r2 in beta))||((r2 in alpha)&&(r1 in beta))) alpha_beta[n]=w
  # coil_glob
  if(((r1 in coil)&&(r2 in coil))) coil_glob[n]=w
  if(((r1 in coil)&&(r2 in alpha))||((r2 in coil)&&(r1 in alpha))) coil_glob[n]=w
  if(((r1 in coil)&&(r2 in beta))||((r2 in coil)&&(r1 in beta))) coil_glob[n]=w
 }
 END{
  # all contacts
  wtot=0.0
  printf("all: COMBINE ARG=C1"); for(i=2;i<=n;i++) printf(",C"i)
  count=0; for(i=1;i<=n;i++) { count++; w=all[i]; wtot+=w; if(count==1) printf(" COEFFICIENTS=%.2f",w); else printf(",%.2f",w) }
  print " PERIODIC=NO # NMAX=",wtot
  # ntail_glob
  wtot=0.0
  count=0; for(i=1;i<=n;i++) { if(i in ntail_glob){ count++; if(count==1) printf("ntail_glob: COMBINE ARG=C"i); else printf(",C"i) } }
  count=0; for(i=1;i<=n;i++) { if(i in ntail_glob){ count++; w=ntail_glob[i]; wtot+=w; if(count==1) printf(" COEFFICIENTS=%.2f",w); else printf(",%.2f",w) } }
  if(count>0) print " PERIODIC=NO # NMAX=",wtot
  # ctail_glob
  wtot=0.0
  count=0; for(i=1;i<=n;i++) { if(i in ctail_glob){ count++; if(count==1) printf("ctail_glob: COMBINE ARG=C"i); else printf(",C"i) } }
  count=0; for(i=1;i<=n;i++) { if(i in ctail_glob){ count++; w=ctail_glob[i]; wtot+=w; if(count==1) printf(" COEFFICIENTS=%.2f",w); else printf(",%.2f",w) } }
  if(count>0) print " PERIODIC=NO # NMAX=",wtot
  # beta_beta
  wtot=0.0
  count=0; for(i=1;i<=n;i++) { if(i in beta_beta){ count++; if(count==1) printf("beta_beta: COMBINE ARG=C"i); else printf(",C"i) } }
  count=0; for(i=1;i<=n;i++) { if(i in beta_beta){ count++; w=beta_beta[i]; wtot+=w; if(count==1) printf(" COEFFICIENTS=%.2f",w); else printf(",%.2f",w) } }
  if(count>0) print " PERIODIC=NO # NMAX=",wtot
  # alpha_alpha
  wtot=0.0
  count=0; for(i=1;i<=n;i++) { if(i in alpha_alpha){ count++; if(count==1) printf("alpha_alpha: COMBINE ARG=C"i); else printf(",C"i) } }
  count=0; for(i=1;i<=n;i++) { if(i in alpha_alpha){ count++; w=alpha_alpha[i]; wtot+=w; if(count==1) printf(" COEFFICIENTS=%.2f",w); else printf(",%.2f",w) } }
  if(count>0) print " PERIODIC=NO # NMAX=",wtot
  # alpha_beta
  wtot=0.0
  count=0; for(i=1;i<=n;i++) { if(i in alpha_beta){ count++; if(count==1) printf("alpha_beta: COMBINE ARG=C"i); else printf(",C"i) } }
  count=0; for(i=1;i<=n;i++) { if(i in alpha_beta){ count++; w=alpha_beta[i]; wtot+=w; if(count==1) printf(" COEFFICIENTS=%.2f",w); else printf(",%.2f",w) } }
  if(count>0) print " PERIODIC=NO # NMAX=",wtot
 }
'

