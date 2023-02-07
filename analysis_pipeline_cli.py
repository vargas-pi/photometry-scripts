"""
This script is a command line interface for running the fiber photometry analysis pipeline
"""


from analysis_pipeline import *
import sys

def input_f(s):
    sys.stdout.flush()
    return input(s)

norm_methods = [norm_to_median_pre_stim, norm_to_405, zscore, zscore_mod]
detrend_methods = [detrend_405, detrend_405_constrained]
def select_norm():
    print('Choose one of the normalization methods:')
    print('1. normalize to the median of pre-stimulus data')
    print('2. normalize to the 405')
    print('3. z-score')
    print('4. z-score to the mean and standard deviation of the pre-stimulus data')
    print('')
    norm_method = input_f('(input the number of the desired choice):')
    return norm_methods[int(norm_method)-1]

def select_detrend():
    print('Choose one of the following detrending methods:')
    print('1. linear fit of the 405 to the 490')
    print('2. linear fit of the 405 to the 490 with positive constraint on the slope')
    print('')
    detrend_method = input_f('(input the number of the desired choice):')
    return detrend_methods[int(detrend_method)-1]

print('==========================================================================')
print('=================fiber-photometry analysis pipeline=======================')
print('============================Alhadeff Lab==================================')
print('==================Monell Chemical Senses Center===========================')
print('==========================================================================')
print(' ')

ans = int(input_f('Would you like to 1. load an existing analysis or 2. start a new analysis? [1/2]: '))

if ans == 1:
    #load in the analysis from the specified file in the output folder
        sys.stdout.flush()
        fname = input_f('Please enter the path to an exported analysis file: ')
        a = load_analysis(fname)
else:
    #ask for necessary parameters for the analysis
    
    t_endrec = input_f('Enter length of recording in seconds from the stimulus time: ')
    t_endrec = float(t_endrec)
    norm_method = select_norm()
    detrend = False
    detrend_method = detrend_405_constrained
    if norm_method != norm_to_405:
        detrend = input_f('Would you like to detrend the data [y/n]: ')
        detrend = detrend not in ['n','no']
        if detrend:
            detrend_method = select_detrend()

    spec_exc_crit=input_f('Would you like to specify the exclusion criteria? (i.e. # of st devs above or below the mean of the data beyond which to exclude) [y/n]: ')
    if spec_exc_crit.lower() in ['n','no']:
        ex=4
    else:
        ex=int(input_f('How many st. devs above or below the mean would you like to define as the limits of the data?: '))
    
    spec_prestim=input_f('Would you like to specify the amount of time pre-stimulus to keep (the default will be 5 minutes)?[y/n]: ')
    if spec_prestim.lower() in ['n','no']:
        t_prestim=300
    else:
        t_prestim=int(input_f('How many seconds pre-stimulus would you like to keep?: '))
    
    a=analysis(norm_method,t_endrec,ex=ex, t_prestim=t_prestim, detrend = detrend, detrend_method=detrend_method)

def load_append_save_cli():

    appending=True

    while appending:
        file=Path(input_f('Please enter the address to an exported data file: ')).resolve()

        if file.suffix=='.json':
            with open(file,'r') as f:
                print('loading data from file...')
                d=json.load(f,object_hook=hook)
        elif file.suffix=='.npy':
            d=np.load(file,allow_pickle=True).tolist()
        else:
            raise Exception('Unrecognized File Format!')

        for i in range(len(d)): #loop through all animal data stored in this run 
            #allow the user to specify the stimulus time in seconds relative to the start of the recording
            if hasattr(d[i],'t_stim'):
                if d[i].t_stim is None:
                    d[i].t_stim=float(input_f(f't_stim for mouse {d[i].mouse_id}:'))
            else:
                d[i].t_stim=float(input_f(f't_stim for mouse {d[i].mouse_id}:'))


            #TODO: finish this up
            if hasattr(d[i],'cond'):
                if d[i].cond is None:
                    resp=input_f(f"if you would like to name the condition for this recording enter it here, otherwise type 'no':  ")
                    if not resp.lower() in ['n','no']:
                        d[i].cond=resp
                    else:
                        d[i].cond=0
            else:
                resp=input_f(f"if you would like to name the condition for this recording enter it here, otherwise type 'no':  ")
                if not resp.lower() in ['n','no']:
                    d[i].cond=resp
                else:
                    d[i].cond=0
                   
            if len(a.raw_data)>0:
                try: d[i].trial=len(a.raw_data.sort_index().loc[d[i].cond,d[i].mouse_id, :])
                except KeyError: d[i].trial=0
            else: d[i].trial=0

            sys.stdout.flush()
            _=a.normalize_downsample(d[i], plot=True)
            #allow the user to decide if they'd like to keep the data and store the data in the appropriate place
            
            try:
                resp=int(input_f('Would you like to 1. keep or 2. discard this data?[1/2]: '))
            except ValueError:
                print ('Answer must be 1 or 2. Try again')
                resp=int(input_f('Would you like to 1. keep or 2. discard this data?[1/2]: '))
            if resp==1:
                a.raw_data[d[i].cond, d[i].mouse_id, d[i].trial] = d[i]
                a.loaded=True
            else:
                a.excluded_raw[d[i].cond, d[i].mouse_id, d[i].trial] = d[i]

        
        cont=input_f('Would you like to continue adding?[y/n]').lower()
        if cont in ['n','no'] :
            appending=False
    a.compute()
    a.save()

def update_params_cli():

    print('Choose one of the following parameters to update:')
    print('1. Length of the recording')
    print('2. Normalization method')
    print('3. Exclusion criteria')
    print('4. Change File Format')
    print('5. Change Pre-Stimulus Time')
    print('6. Whether or not to detrend')
    print('7. Set the detrend method')
    print('8. Set the downsampling frequency')
    print('')
    choice=input_f('(input the number of the desired choice):')

    if choice=='1':
        a.t_endrec = int(input_f( 'Enter length of recording in seconds from the stimulus time:' ))
    elif choice=='2':
        a.norm_method = select_norm()
    elif choice=='3':
        a.ex=int(input_f( 'How many st. devs above or below the mean would you like to define as the limits of the data?: ' ))
    elif choice=='4':
        a.file_format=['npy','json'][-1+int(input_f( 'What file format would you like? 1. npy, 2. json [1/2]: ' ))]
    elif choice=='5':
        a.t_prestim=int(input_f( 'Enter the desired pre-stimulus time in seconds:' ))
    elif choice=='6':
        detrend = input_f('Would you like to detrend the data [y/n]: ')
        a.detrend = detrend not in ['n','no']
    elif choice=='7':
        a.detrend_method = select_detrend()
    elif choice=='8':
        a.ds_freq = int(input_f( 'What frequency in Hz would you like the data downsampled to?:' ))
    a.compute()

def remove_mouse_cli():
    a.remove_mouse(input_f("Enter the id of the mouse you'd like to remove: "))

def retrieve_excluded_cli():
    m = input_f("Enter the id of the mouse you'd like to retrieve: ")
    c = input_f("please enter the condition you would like to retrieve from this mouse's data (to retrieve all type 'all'): ")
    if c=='all':
        c=None
    a.retrieve_excluded(m, c)

def plot_both_cli():
    ans = input_f('Would you like to specify the y limit of the plot? [y/n]: ')
    if ans.lower() not in ['n','no']:
        l = float(input_f('Please specify the lower bound of the plot: '))
        u = float(input_f('Please specify the upper bound of the plot: '))
        ylim = (l,u)
    else:
        ylim = None
    a.plot_both(ylim = ylim)

def plot_490_cli():
    ans = input_f('Would you like to specify the y limit of the plot? [y/n]: ')
    if ans.lower() not in ['n','no']:
        l = float(input_f('Please specify the lower bound of the plot: '))
        u = float(input_f('Please specify the upper bound of the plot: '))
        ylim = (l,u)
    else:
        ylim = None
    a.plot_490(ylim = ylim)
        

def bin_plot_cli():
    binsize = int(input_f('How big, in seconds, would you like the bins: '))
    ans = input_f('Would you like to specify the y limit of the plot? [y/n]: ')
    if ans.lower() not in ['n','no']:
        l = float(input_f('Please specify the lower bound of the plot: '))
        u = float(input_f('Please specify the upper bound of the plot: '))
        ylim = (l,u)
    else:
        ylim = None
    a.bin_plot(binsize, save=True, ylim = ylim)

def bin_auc_cli():
    start=int(input_f('Enter the beginning of the period in seconds relative to the stimulus onset: '))
    end=int(input_f('Enter the end of the period: '))
    a.bin_auc(start,end)

def bin_avg_cli():
    start=int(input_f('Enter the beginning of the period to average in seconds relative to the stimulus onset: '))
    end=int(input_f('Enter the end of the period: '))
    a.bin_avg(start,end,save=True)

def ind_peak_df_f_cli():
    ans=int(input_f('Would you like to compute the 1. max or 2. min? [1/2] '))-1
    opts=['max','min']
    vals=int(input_f('Would you like the 1. raw values or 2. integral values? '))-1
    methods=['False', 'True']
    a.ind_peak_df_f(opts[ans],save=True, trap=methods[vals])

def mean_peak_df_f_cli():
    ans=int(input_f('Would you like to compute the 1. max or 2. min? [1/2] '))-1
    opts=['max','min']
    vals=int(input_f('Would you like the 1. raw values or 2. integral values? '))-1
    methods=['False', 'True']
    a.mean_peak_df_f(opts[ans],save=True, trap=methods[vals])

def time_to_half_pk_cli():
    ans=int(input_f('Would you like to compute the 1. max or 2. min? [1/2] '))-1
    opts=['max','min']
    a.time_to_half_pk(opts[ans], filtered=True, pr=True, save=True)

def plot_ind_trace_cli():
    m=input_f('please enter the name of the mouse you would like to view: ')
    c=input_f("please enter the condition you would like to plot this mouse's data from (to plot all type 'no'): ")
    c=None if c.lower() in ['n','no'] else c
    p4=input_f('would you like to plot the 405? [y/n]: ').lower() in ['y','ye','yes']
    a.plot_ind_trace(m,cond=c,plot_405=p4)

running=True
if not a.loaded:
    load_append_save_cli()
    a.plot_both()

while running:
    print('------------------------------------')
    print('Choose one of the following tasks:')
    print('1. plot 490 and 405')
    print('2. plot just the 490')
    print('3. add data to this analysis')
    print('4. save this analysis')
    print('5. find the peak ∆f/f (min/max) for individual mice')
    print('6. find the peak ∆f/f (min/max) for individual mice at the location of the peak in the mean signal')
    print('7. find the time to half peak')
    print('8. find the average value over a specified portion of data')
    print('9. find the area under the curve for a specified portion of data')
    print('10. bin and plot the data')
    print('11. remove a mouse from this analysis')
    print('12. retrieve an excluded mouse from this analysis')
    print('13. update the parameters of this analysis')
    print('14. export normalized 490 data to .mat')
    print('15. plot the trace for an individual mouse')
    print('16. exit')
    print('')
    ans=input_f('(input the number of the desired task): ')



    tasks={
        '1':plot_both_cli,
        '2':plot_490_cli,
        '3':load_append_save_cli,
        '4':a.save,
        '5':ind_peak_df_f_cli,
        '6':mean_peak_df_f_cli,
        '7':time_to_half_pk_cli,
        '8':bin_avg_cli,
        '9':bin_auc_cli,
        '10':bin_plot_cli,
        '11':remove_mouse_cli,
        '12':retrieve_excluded_cli,
        '13':update_params_cli,
        '14':a.export_to_mat,
        '15':plot_ind_trace_cli
        }
    
    try:
        tasks[ans]()
    except KeyError:
        running=False
