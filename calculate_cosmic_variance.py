#------------------------------------------------Import required packages---------------------------------------
import numpy
from multiprocessing import Pool
import pickle
import scipy.integrate
import scipy.interpolate
#-----------------------------Internal variables: DO NOT CHANGE---------------------------------------------
arc_sec_to_radians=1./(180/3.14)*(1./3600.)
#-------------------------------Read fiducial cosmology---------------------------------------------
f=open('cosmology.txt')
lines=f.readlines()
for line in lines:
    if ('-------' not in line):
        #line= line.replace('=',' ')
        #line=line.split()
        #print line
        exec(line)
        
#-------------------------------Read default inputs for the monte carlo integration---------------------------------------------

f=open('parameters_for_monte_carlo.txt')
lines=f.readlines()
for line in lines:
    if ('-------' not in line):
            #line= line.replace('=',' ')
            #line=line.split()
            #print line
            exec(line)



#---------------------------Function definitions-----------------------------------------
H = lambda z: H0*(om0*(1+z)**(3) + oml)**0.5

def box_size_to_delta_z(z_mid,BOX_LENGTH): # Calculates the box_size in Mpc/h for a given redshift interval 
    def er(a,zz):
        return numpy.subtract(BOX_LENGTH,scipy.integrate.quad(lambda z:c/H(z),zz-a/2,zz+a/2))[0]
    return scipy.optimize.newton(er,2,args=(z_mid,))

def delta_z_to_box_size(z_mid,delta_z): # Calculates the redshift interval for a given box_size in Mpc/h 
    box_size=scipy.integrate.quad(lambda z:c/H(z),z_mid-delta_z/2,z_mid+delta_z/2)[0]
    return box_size

def get_survey_dimensions(central_redshift,redshift_width,survey_length_in_arc_sec): # Calculates the survey dimensions in cMpc/h for a given field of view (in arc-sec) and redshift width 
    distance_to_central_redshift=comoving_distance(0,central_redshift)
    #print "comoving_distance:",distance_to_central_redshift/h
    survey_length_in_radians=survey_length_in_arc_sec*arc_sec_to_radians
    survey_length_in_cmpc_h=survey_length_in_radians*distance_to_central_redshift
    redshift_width_in_Mpc=delta_z_to_box_size(central_redshift,redshift_width)
    survey_dimensions=numpy.array([survey_length_in_cmpc_h,survey_length_in_cmpc_h*survey_aspect_ratio,redshift_width_in_Mpc])
    return survey_dimensions
    

def comoving_distance(z_low,z_up): # Calculates the comoving distance between two redshifts 
    box_size=scipy.integrate.quad(lambda z:c/H(z),z_low,z_up)[0]
    return box_size

#-----------------------------------------------------------------------------------------------------------------------
if Integrate_over_logspace:
    def test_function(input_array):
        modified_input_array=10**input_array
        xi=numpy.product(modified_input_array)**N
        #xi=3.
        return xi*numpy.product(modified_input_array)
    
    def correlation_function(input_array):
        modified_input_array=10**input_array
        x1,y1,z1,x2,y2,z2=modified_input_array
        r=numpy.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
        #print numpy.sqrt(3.)*x2,numpy.sqrt(3.)*y2,numpy.sqrt(3.)*z2,r
        xi=10**power_law(numpy.log10(r),fitted_exponent,fitted_normalization)
        return xi*numpy.product(modified_input_array)
    
    def linear_dark_matter(input_array):
        modified_input_array=10**input_array
        x1,y1,z1,x2,y2,z2=modified_input_array
        r=numpy.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
        #print numpy.sqrt(3.)*x2,numpy.sqrt(3.)*y2,numpy.sqrt(3.)*z2,r
        xi=10**gen_dm_corr(numpy.log10(r))
        if ((NON_LINEAR_POWER_SPECTRUM==1)&(r<10.)):
            xi=10**gen_dm_corr_nl(numpy.log10(r))
        #xi=3.    
        return xi*numpy.product(modified_input_array)
    
else:
    def test_function(input_array):
        return numpy.product(input_array)**N
    
    def correlation_function(input_array):
        modified_input_array=input_array
        x1,y1,z1,x2,y2,z2=modified_input_array
        r=numpy.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
        #print numpy.sqrt(3.)*x2,numpy.sqrt(3.)*y2,numpy.sqrt(3.)*z2,r
        xi=10**power_law(numpy.log10(r),fitted_exponent,fitted_normalization)
        return xi
    
    def linear_dark_matter(input_array):
        modified_input_array=10**input_array
        x1,y1,z1,x2,y2,z2=modified_input_array
        r=numpy.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
        #print numpy.sqrt(3.)*x2,numpy.sqrt(3.)*y2,numpy.sqrt(3.)*z2,r
        xi=10**gen_dm_corr(numpy.log10(r))
        if ((NON_LINEAR_POWER_SPECTRUM==1)&(r<10.)):
            xi=10**gen_dm_corr_nl(numpy.log10(r))
        #xi=3.
        return xi

def power_law(log_distance,exponent,normalization):
    return exponent*log_distance+normalization
    
def get_test_integral():
    #print "Check Final values:",Final_values
    #print "Check Initial values:",Initial_values
    if (Integrate_over_logspace):
        Initial_values_temp=numpy.log10(Initial_values)
        Final_values_temp=numpy.log10(Final_values)
    else:
        Initial_values_temp=Initial_values
        Final_values_temp=Final_values
    #print Initial_values_temp
    #print Final_values_temp
    print numpy.product(((10**Final_values_temp)**(N+1)-(10**Initial_values_temp)**(N+1))/(N+1))
    widths=Final_values_temp-Initial_values_temp
    scaled_points=numpy.random.rand(Number_of_points,Dimensions)
    actual_points=scaled_points*widths+Initial_values_temp
    #---------------------------------Testing the range----------------------------------
    print "Testing the range"
    tolerance=1e-2
    t_actual_points=numpy.transpose(actual_points)
    for i in range(0,Dimensions):
        left_edge_gap=numpy.abs(numpy.amin(t_actual_points[i])-Initial_values_temp[i])
        right_edge_gap=numpy.abs(numpy.amax(t_actual_points[i])-Final_values_temp[i])
        if (left_edge_gap>tolerance):
            print "Warning! Range may not be correct for axis ",i," with left_edge_gap ",left_edge_gap
        if (right_edge_gap>tolerance):
            print "Warning! Range may not be correct for axis ",i," with right_edge_gap ",right_edge_gap
    #--------------------------------Performing a test integration-----------------------------
    if (TEST_MODE):
        if (PARALLEL_PROCESSORS):
            print "parallel calculation"
            p=Pool(Number_of_processors)
            values=numpy.array(p.map(test_function,actual_points)) 
            p.close()
        else:
            print "serial calculation"
            
            values=numpy.array([test_function(point) for point in actual_points]) 
    
        phase_space_volume=numpy.product(widths)
        integral=numpy.average(values)*phase_space_volume
    
    
    
        if (Integrate_over_logspace):
            integral*=numpy.log(10)**Dimensions
            analytical_integral=numpy.product(((10**Final_values_temp)**(N+1)-(10**Initial_values_temp)**(N+1))/(N+1))
        else:
            analytical_integral=numpy.product((Final_values_temp**(N+1)-Initial_values_temp**(N+1))/(N+1))
        return integral,analytical_integral
    
    else:
        if (PARALLEL_PROCESSORS):
            print "Performing parallel calculation"
            p=Pool(Number_of_processors)
            
            if (LINEAR_DARK_MATTER_MODE):
                values=numpy.array(p.map(linear_dark_matter,actual_points)) 
                p.close()   
            else:    
                values=numpy.array(p.map(correlation_function,actual_points)) 
                p.close()
        else:
            print "Performing serial calculation"
            values=numpy.array([correlation_function(point) for point in actual_points]) 

        #print "Volume averaged correlation function", numpy.average(values)
    
        phase_space_volume=numpy.product(widths)
        integral=numpy.average(values)*phase_space_volume
        if (Integrate_over_logspace):
            integral*=numpy.log(10)**Dimensions
        return integral
    
#-----------------------------OVERRIDING DEFAULT INPUTS------------------------------------------------
TEST_MODE=False             #If this is True, it calculates cosmic variance for an analytically integrable function 
PARALLEL_PROCESSORS=True
LINEAR_DARK_MATTER_MODE=True   #IF this is true, it calculates the cosmic variance for dark matter-----
NON_LINEAR_POWER_SPECTRUM=0
#-----------------------------Reading the survey input and the associated parameters-------------------

survey_type='COSMOS'
prefix='survey_specifications_'
f=open(prefix+survey_type+'.txt')
lines=f.readlines()
for line in lines:
    if ('-------' not in line):
            #print line
            exec(line)  
            
            
Number_of_points_space=[10**5,10**6,10**7][:]
for Number_of_points in [10**8]:
        
        integral_space=[]
        redshift_for_plot=[]
        N=-2        
        for redshift in [0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5]:
        #for redshift in [0.7]:
            #try:
                print redshift
                r_dm,xi_dm=pickle.load(open('../project_3/linear_matter_correlation/z%.2f_correlation_second_round.pickle'%(redshift)))                
                mask_dm=xi_dm>0
                gen_dm_corr=scipy.interpolate.interp1d(numpy.log10(r_dm[mask_dm]),numpy.log10(xi_dm[mask_dm]),fill_value='extrapolate')                
                dimensions=get_survey_dimensions(redshift,redshift_width,survey_length_in_arc_sec)
                Final_values=numpy.append(dimensions,dimensions)
                smallest_scale=get_survey_dimensions(redshift,redshift_width,survey_resolution_in_arc_sec)[0]
                Dimensions=len(Final_values)
                Initial_values=numpy.array([smallest_scale]*Dimensions)
                survey_volume=numpy.product(dimensions)
                integrals=get_test_integral()
                Phase_space_volume_sq=numpy.product(Final_values-Initial_values)                
                cosmic_variance=numpy.sqrt(integrals/survey_volume**2)
                if (LINEAR_DARK_MATTER_MODE):
                    pickle.dump([cosmic_variance],open('./cosmic_variance_outputs/%s_redshift_%.2f_log_no_of_points_%d'%(survey_type,redshift,numpy.log10(Number_of_points)),'w'))
                print cosmic_variance
            #except:
            #    print "No file"
            #    aaa=1







