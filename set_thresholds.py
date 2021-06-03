
#!/opt/annaconda/python3
from scipy.stats import binom
import sys



if len(sys.argv)>1:
  site_freq_file=open(sys.argv[1],'r')
else:
  try:
    site_freq_file=open('site_frequency.txt','r')
  except:
    sys.exit('No site frequency file specified and default file cannot be found')
if len(sys.argv)>2:
  f=open(sys.argv[2],'r')
else:
  try:
    f=open('threshold_parameters.txt','r')
  except:
    sys.exit('No parameter file specified and default file cannot be found')
line=f.readline()
while line:
  line_list=line.split('=')
  if line_list[0]=='external_haploid_sample_size':
    sample_size1=int(line_list[1].strip('\n'))
  elif line_list[0]=='internal_haploid_sample_size':
    sample_size2=int(line_list[1].strip('\n'))
  elif line_list[0]=='sensitivity_thresholds':
    sensitivity_list=[float(item) for item in line_list[1].strip('\n').split(',')]
  elif line_list[0]=='internal_to_external_relative_risk':
    rr=float(line_list[1])
  elif line_list[0]=='max_internal_allele_frequency':
    max_internal_allele_frequency=float(line_list[1].strip('\n'))
  line=f.readline()
try:
  check=sample_size1+sample_size2+sum(sensitivity_list)+rr+max_internal_allele_frequency
except:
  sys.exit('One or more required parameters are undefined')
f.close()
if len(sys.argv)>3:
  output_file=sys.argv[3]
else:
  output_file='or_filters.txt'
try:
  g=open(output_file,'w')
except:
  sys.exit('Cannot create file '+output_file)
sum=0.0
k={}
line=site_freq_file.readline()
line=site_freq_file.readline()
while line:
  line_list=line.split()
  k[int(line_list[0].strip('"'))]=int(line_list[1].strip('"'))
  line=site_freq_file.readline()
site_freq_file.close()
for i,count in k.iteritems():
  sum+=i
for i in xrange(sample_size1):
  if i in k:
    k[i]=k[i]/float(sum)
  else:
    k[i]=0.0
g.write("site_count\tmax_accuracy\tor_to_maximize_accuracy")
for s in sensitivity_list:
  g.write('\t'+str(s))
for i in xrange(1,int(sample_size2*max_internal_allele_frequency)):
  #g.write('\n'+str(i))
  line=''
  null_site_freq_cond_on_x={}
  alt_site_freq_cond_on_x={}
  null_prop_sum_current=0.0
  alt_prop_sum_current=0.0
  null_prop_sum=0.0
  alt_prop_sum=0.0
  null_site_freq_cond_on_x[i]={}
  alt_site_freq_cond_on_x[i]={}
  for j in xrange(0,sample_size1-sample_size2): 
    null_p=float(i+j)/float(sample_size1)
    term=rr*sample_size2/(sample_size1-sample_size2)
    alt_p=term*(i+j)/(sample_size2*(1+term))
    null_site_freq_cond_on_x[i][j]=binom.pmf(i,sample_size2,null_p)*k[i+j] 
    alt_site_freq_cond_on_x[i][j]=binom.pmf(i,sample_size2,alt_p)*k[i+j] 
  null_prop_sum=0.0
  for j,prop in null_site_freq_cond_on_x[i].iteritems():
    if prop>0:
      null_prop_sum+=prop 
  alt_prop_sum=0.0
  for j,prop in alt_site_freq_cond_on_x[i].iteritems():
    if prop>0:
      alt_prop_sum+=prop 
  min_odds_ratio=[0]*len(sensitivity_list)
  min_odds_ratio_acc=0
  max_accuracy=0
  sensitivity_index=0
  for j,null_prop in null_site_freq_cond_on_x[i].iteritems():
    alt_prop=alt_site_freq_cond_on_x[i][j]
    if null_prop>0:
      if j>0:
        odr=(float(i)/sample_size2)/(j/float((sample_size1-sample_size2)))
      else:
        odr=(float(i)/sample_size2)/(0.5/float((sample_size1-sample_size2)))
      null_prop_sum_current+=null_prop/null_prop_sum
      alt_prop_sum_current+=alt_prop/alt_prop_sum
      while sensitivity_index<len(sensitivity_list) and alt_prop_sum_current>=sensitivity_list[sensitivity_index] and min_odds_ratio[sensitivity_index]==0:
        #g.write('\t'+str(odr))
        line+='\t'+str(odr)
        min_odds_ratio[sensitivity_index]=odr
        sensitivity_index+=1
      if (1-null_prop_sum_current)+alt_prop_sum_current>max_accuracy:
        max_accuracy=(1-null_prop_sum_current)+alt_prop_sum_current
        min_odds_ratio_acc=odr
  g.write('\n'+str(i)+'\t'+str(max_accuracy/2.0)+'\t'+str(min_odds_ratio_acc)+line) 
  print str(round(100.0*i/int(sample_size2*max_internal_allele_frequency),2))+'% complete'
g.close()
