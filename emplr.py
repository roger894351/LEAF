from scipy.stats import binom
import math

def init(site_freq_file):
  global sum,k
  sum=0.0
  k={}
  f=open(site_freq_file,'r')
  line=f.readline()
  line=f.readline()
  while line:
    line_list=line.split()
    k[int(line_list[0].strip('"'))]=int(line_list[1].strip('"'))
    line=f.readline()
  for i,count in k.iteritems():
    sum+=i
  for i in k.keys():
    k[i]=k[i]/float(sum)

def get_likelihood(i,j,rr,internal_sample_size,external_sample_size,significant_digits=2):
  term=float(rr)*internal_sample_size/(external_sample_size-internal_sample_size)
  p=term*(i+j)/(internal_sample_size*(1+term))
  if i+j not in k:
    k[i+j]=0
  l=binom.pmf(i,internal_sample_size,p)*k[i+j] 
  if l==0:
    return l
  else:
    return round(l, significant_digits - int(math.floor(math.log10(abs(l)))) - 1)

def find_ml(i,j,internal_sample_size,external_sample_size,significant_digits=2):
  if j==0:
    mid_rr=(float(i)/internal_sample_size)/(0.5/external_sample_size)
  else:
    mid_rr=(float(i)/internal_sample_size)/(float(j)/external_sample_size)
  max_rr=mid_rr*100
  min_rr=mid_rr/100
  max_l=get_likelihood(i,j,max_rr,internal_sample_size,external_sample_size,significant_digits)
  min_l=get_likelihood(i,j,min_rr,internal_sample_size,external_sample_size,significant_digits)
  mid_l=min_l
  while min_l!=max_l:
    mid_l=get_likelihood(i,j,mid_rr,internal_sample_size,external_sample_size,significant_digits) 
    print min_l,min_rr,max_l,max_rr,mid_l,mid_rr
    lower_mid=(mid_rr+min_rr)/2.0
    upper_mid=(mid_rr+max_rr)/2.0
    lower_mid_l=get_likelihood(i,j,lower_mid,internal_sample_size,external_sample_size,significant_digits)
    upper_mid_l=get_likelihood(i,j,upper_mid,internal_sample_size,external_sample_size,significant_digits)
    while upper_mid_l==lower_mid_l and upper_mid_l==mid_l and not ((upper_mid_l==max_l and max_l>=min_l) or (lower_mid_l==min_l and min_l>=max_l)):
      lower_mid=(lower_mid+min_rr)/2.0
      upper_mid=(upper_mid+max_rr)/2.0
      lower_mid_l=get_likelihood(i,j,lower_mid,internal_sample_size,external_sample_size,significant_digits)
      upper_mid_l=get_likelihood(i,j,upper_mid,internal_sample_size,external_sample_size,significant_digits)
    if mid_l>=upper_mid_l and mid_l>=lower_mid_l:
      max_rr=upper_mid
      min_rr=lower_mid
      max_l=upper_mid_l
      min_l=lower_mid_l
    elif upper_mid_l>lower_mid_l:
      min_rr=lower_mid
      min_l=lower_mid_l
      mid_rr=(min_rr+max_rr)/2.0
      mid_l=get_likelihood(i,j,mid_rr,internal_sample_size,external_sample_size,significant_digits)
    elif lower_mid_l>upper_mid_l:
      max_rr=upper_mid
      max_l=upper_mid_l
      mid_rr=(min_rr+max_rr)/2.0
      mid_l=get_likelihood(i,j,mid_rr,internal_sample_size,external_sample_size,significant_digits)
    else:
      sys.exit('Unexpected result in likelihood search, upper_mid_rr='+str(upper_mid_rr)+' lower_mid_rr='+str(lower_mid_rr)+' mid_rr='+str(mid_rr))
  return mid_rr,mid_l

def get_lr(i,j,internal_sample_size,external_sample_size,significant_digits=2):
  null_l=get_likelihood(i,j,1.0,internal_sample_size,external_sample_size,significant_digits)
  rr,alt_l=find_ml(i,j,internal_sample_size,external_sample_size,significant_digits)
  if rr>1:
    lr=alt_l/null_l
  else:
    lr=null_l/alt_l
  return lr
