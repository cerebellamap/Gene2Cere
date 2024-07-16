import os
import pandas as pd
import numpy as np
import nibabel as nib
from glob import glob
from joblib import Parallel, delayed
from sklearn.model_selection import KFold
from sklearn.cross_decomposition import PLSRegression
from sklearn import model_selection, linear_model
from sklearn.metrics import mean_squared_error
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import matplotlib as mpl
from tqdm import tqdm
import time
start = time.time()
import multiprocessing
import subprocess 
import sys
sys.path.append('/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Toolbox/')
from abagen_cere import allen
from abagen_cere import datasets

# colorbar
RdBu_200 = cm.get_cmap('RdBu', 200)(np.arange(200))
cm.get_cmap('RdBu', 200)
RdBu_Hex = [mpl.colors.to_hex(x) for x in np.flipud(cm.get_cmap('RdBu')(np.linspace(0,1,200)))[:,0:3]]
RdBu_Hex = RdBu_Hex[: :-1]   

from PIL import Image, ImageDraw, ImageFont
def get_dominant_colors(infile):
    image = Image.open(infile)
    
    # 缩小图片，否则计算机压力太大
    small_image = image.resize((80, 80))
    result = small_image.convert(
        "P", palette=Image.ADAPTIVE, colors=10
    )  
	
	# 10个主要颜色的图像

    # 找到主要的颜色
    palette = result.getpalette()
    color_counts = sorted(result.getcolors(), reverse=True)
    colors = list()

    for i in range(10):
        palette_index = color_counts[i][1]
        dominant_color = palette[palette_index * 3 : palette_index * 3 + 3]
        colors.append(tuple(dominant_color))

    # print(colors)
    return colors

def RGB_to_Hex(rgb):
    RGB = rgb.split(',')
    color = '#'
    for i in RGB:
        num = int(i)
        color += str(hex(num))[-2:].replace('x', '0').upper()
    return color


image_path = "/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/colorbar/Hex_1.jpg"
color_1 = get_dominant_colors(image_path)
# print(color_1)
color_rgb = [str(x).replace('(','').replace(')','')  for x in color_1]
Hex_1 = [RGB_to_Hex(x) for x in color_rgb]

image_path = "/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/colorbar/Hex_2.jpg"
color_2 = get_dominant_colors(image_path)
# print(color_2)
color_rgb = [str(x).replace('(','').replace(')','')  for x in color_2]
Hex_2 = [RGB_to_Hex(x) for x in color_rgb]


jet_10 = cm.get_cmap('jet', 10)(np.arange(10))
jet_9 = cm.get_cmap('jet', 9)(np.arange(9))
RdBu_9 = cm.get_cmap('RdBu', 9)(np.arange(9))



# settings.py
class Settings:
    def __init__(self, data_dir=None, output_dir=None):
        self.data_dir = data_dir or '/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Data/'
        self.output_dir = output_dir or '/n02dat01/users/ypwang/Gradient/STARProtocol/Try/Output/'
# Initialize settings
global_settings = Settings()


def Step01_Input(input_file_name, data_dir=None, output_dir=None):
    # Directories
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    input_file = os.path.join(data_dir, input_file_name)    
    ahba_dir = f'{data_dir}AHBA_data/'  
    # Check for the existence of the microarray folder and fetch data if not present
    if not os.path.exists(os.path.join(ahba_dir, 'microarray')):
        print('Fetching microarray data...')
        donors = ['12876', '14380', '15496', '15697', '9861', '10021']  # Example donor IDs
        files = datasets.fetch_microarray(data_dir=ahba_dir, donors=donors, verbose=1, n_proc=4)
  
                               
    coord_pth = f'{data_dir}AHBA_data/Cerebellum_sample_coordinates_corrected/'

    # Supporting functions
    def disbtw(surf, aa):
        surf = np.array(surf)
        aa = np.expand_dims(aa, axis=0)
        dis = np.sqrt(((surf - aa) ** 2).sum(axis=1))
        return dis

    def Grid2world(GridList, affine):
        GridList = np.array(GridList)
        affine = np.array(affine)
        Grid_coord = np.concatenate([GridList, np.ones([len(GridList), 1])], axis=1)
        MNI_coord = np.dot(Grid_coord, affine.T).round(2)
        return MNI_coord[:, :3]


    # Get the cerebellar sample-wise transcriptom matrix
    print('Get the cerebellar sample-wise transcriptom matrix')
    expression_all, report_all = allen.get_samples_in_mask_annot(data_dir=ahba_dir, corrected_mni=True, return_report=True, Allsample_nodrop='noB')
    report = pd.DataFrame(report_all[report_all.slab_type=='CB'])
    report = report[~report.structure_name.str.contains('nucleus')]
    expression = expression_all[expression_all.index.isin(report.index)]

    # Read in the corrected coordinates
    print('Read in the corrected coordinates')
    sasheets = sorted(glob(os.path.join(coord_pth,'*/*_mni_coordinates.csv')))
    ref = []
    for sheet in sasheets:
        did = sheet.split('/')[-2]
        sa = pd.read_csv(sheet, header=None)
        sa.columns = ['suit_x', 'suit_y', 'suit_z']
        sa.index = report[report.donor==did].index
        ref.append(sa)
    SA = pd.concat(ref)
    report = pd.merge(report, SA, left_index=True, right_index=True)

    # Assign the imaging features to each cerebellar sample
    print('Assign the imaging features to each cerebellar sample')
    Input_file = nib.load(input_file)
    Input = Input_file.get_fdata()
    affine = Input_file.affine
    grid_list = np.nonzero(Input)
    grad_value = Input[grid_list]
    MNI_loc = Grid2world(np.array(grid_list).T, affine)
    report['Y'] = [np.mean(grad_value[disbtw(MNI_loc, loc) < 4]) for loc in report[['suit_x', 'suit_y', 'suit_z']].values.astype(float)]

    # Delete the sample whose imaging features is nan
    report_cur = report[~np.isnan(report['Y'])]
    expression_cur = expression[~np.isnan(report['Y'])]
    report = report_cur.iloc[np.where(report_cur['Y']!=0)]
    expression = expression_cur.iloc[np.where(report_cur['Y']!=0)]

    # Save
    print('Save Gene_expression.csv and Sample_info.csv')
    pd.DataFrame(expression).to_csv(f'{output_dir}Step01_Gene_expression.csv')
    pd.DataFrame(report).to_csv(f'{output_dir}Step01_Sample_info.csv')



def Step02_Comp_eval_run(num_repeat, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    all_best_comps, all_r, all_preds = [], [], []
    for i in range(num_repeat):
        print(f'Repeat {i + 1}/{num_repeat}')
        best_comps, r, preds = Component_eval(i)
        all_best_comps.append(best_comps)
        all_r.append(r)
        all_preds.append(preds)
    np.save(f'{output_dir}/Step02_Comp_eval_run_{num_repeat}x10cv_all_best_comps', all_best_comps)
    np.save(f'{output_dir}/Step02_Comp_eval_run_{num_repeat}x10cv_all_r', all_r) 
    return f'{output_dir}/Step02_Comp_eval_run_{num_repeat}x10cv_all_best_comps', f'{output_dir}/Step02_Comp_eval_run_{num_repeat}x10cv_all_r'

def Step02_Comp_eval_visualization(all_best_comps, all_r, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    all_best_comps = np.load(f'{all_best_comps}.npy', allow_pickle=True)
    all_r = np.load(f'{all_r}.npy', allow_pickle=True)
    df = pd.DataFrame()
    for x in range(len(all_best_comps)):
        print(x)
        df_cur = pd.DataFrame(columns = ["Component_9fold_bymse", "Correlation"])
        df_cur.Component_9fold_bymse = all_best_comps[x]
        df_cur.Correlation = np.array(all_r[x])
        df = pd.concat([df, df_cur], axis = 0)
    r_nc_plot = pd.DataFrame({'cn':df.Component_9fold_bymse, '100*10cv_score':abs(df.Correlation.values)})
    for cn in r_nc_plot.cn.unique():
        print(cn)
        print(np.mean(r_nc_plot.loc[r_nc_plot.cn == cn,'100*10cv_score']))
        print(np.median(r_nc_plot.loc[r_nc_plot.cn == cn,'100*10cv_score']))
    mean_scores = r_nc_plot.groupby('cn')['100*10cv_score'].mean()
    max_mean_cn = mean_scores.idxmax()
    fig = plt.figure(figsize=(12,10))

    my_pal = {cn: RdBu_200[30] if cn == max_mean_cn else RdBu_200[185] for cn in r_nc_plot.cn.unique()}
    # my_pal = {cn: jet_9[7] if cn == 6 else jet_9[1] for cn in r_nc_plot.cn.unique()}
    sns.boxplot(x="cn", y='100*10cv_score', data=r_nc_plot, showfliers=False, palette=my_pal, width = 0.8, showmeans=True, 
                boxprops = { 'edgecolor':'grey'},
                meanprops={"marker":"o", "markersize":"10", "markerfacecolor":"white",  "linewidth":1, "markeredgecolor":"black"},
                medianprops={"color":"grey", "linewidth":2},
                whiskerprops={"color":"grey", "linewidth":2},
                capprops={"color":"grey", "linewidth":2}) # palette="hls"
    sns.swarmplot(x="cn", y='100*10cv_score', data=r_nc_plot, color ='gray',size = 4, alpha=0.66)
    sns.despine()
    plt.ylim(0, 0.8)
    plt.yticks([0,0.2,0.4,0.6,0.8], fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel('Component number', fontsize=25)
    plt.ylabel("Correlation", fontsize=25)

    fig_name = 'Step02_Comp_eval_visualization.png'
    plt.savefig(os.path.join(output_dir, fig_name), dpi=500, bbox_inches='tight')
    plt.show()

def Component_eval_embededcv(tr_ind, te_ind, expression, report, ncomps, random_seed):
    data_dir = global_settings.data_dir
    cv_inner = 10
    sel_inner = KFold(n_splits=cv_inner, shuffle=True, random_state=(123 * (random_seed+1)))
    expression_cur = expression.iloc[tr_ind,:]
    mse = []

    for tr_i, te_i in sel_inner.split(expression_cur):
        tr_cols = expression.index[tr_ind][tr_i]
        tr_set = expression.loc[tr_cols,:]
        tr_y = report['Y'].loc[tr_cols]
        te_cols = expression.index[tr_ind][te_i]
        te_set = expression.loc[te_cols,:]
        te_y = report['Y'].loc[te_cols]

        mse_iter = [mean_squared_error(te_y, PLSRegression(n_components=nc).fit(tr_set, tr_y).predict(te_set).flatten()) for nc in ncomps]
        mse.append(mse_iter)

    mse = np.sum(mse, axis=0)
    nc = ncomps[np.argmin(mse)]
    clf = PLSRegression(n_components=nc)
    tr_set = expression.iloc[tr_ind,:]
    tr_y = report['Y'].iloc[tr_ind]
    te_set = expression.iloc[te_ind,:]
    te_y = report['Y'].iloc[te_ind]
    mod = clf.fit(tr_set, tr_y)
    preds = pd.Series(mod.predict(te_set).flatten(), index=expression.index[te_ind])

    r = stats.pearsonr(mod.predict(te_set).flatten(), te_y)[0]
    return nc, r, preds

def Component_eval(random_seed, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    expression = pd.read_csv(os.path.join(output_dir, 'Step01_Gene_expression.csv'), index_col=0)
    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)
    ncomps = np.arange(1, 11)
    sel = KFold(n_splits=10, shuffle=True, random_state=(123 * (random_seed + 1)))

    results = Parallel(n_jobs=-1)(delayed(Component_eval_embededcv)(tr_ind, te_ind, expression, report, ncomps, random_seed) for tr_ind, te_ind in sel.split(expression))

    best_comps, r, preds = zip(*results)
    preds = pd.concat(preds, axis=1)
    return best_comps, r, preds


def Step02_Model_eval(n_components,cv_repeat, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    expression = pd.read_csv(os.path.join(output_dir, 'Step01_Gene_expression.csv'), index_col=0)
    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)

    clf = PLSRegression(n_components)
    cv_strategy = 10
    
    # Running the parallel processing
    results = Parallel(n_jobs=-1)(delayed(perform_CV)(i, expression, report, clf, cv_strategy, cv_repeat) for i in range(cv_repeat))
    # Extract predictions  from results
    preds_rep = pd.DataFrame([result for result in results]).T
    preds_rep.columns = range(cv_repeat)
    preds_rep.index = expression.index

    score_true = []
    for i in range(cv_repeat):
        score_true.append(stats.pearsonr(preds_rep.iloc[:,i].values, report['Y'])[0])
    score_df = pd.DataFrame(score_true)
    median_index = score_df.sort_values(by=0).index[50] 
    median_score = score_df.sort_values(by=0).iloc[50]

    # Saving results to CSV
    preds_rep.to_csv(os.path.join(output_dir, f'Step02_PLSR_{cv_repeat}x10cv_preds.csv'))
    print(f'median_index: {median_index}, median_score: {median_score[0]}')
    return median_index, median_score[0]

def Step02_Model_eval_visulization(n_components,cv_repeat,median_index,median_score, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    expression = pd.read_csv(os.path.join(output_dir, 'Step01_Gene_expression.csv'), index_col=0)
    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)

    preds_rep = pd.read_csv(os.path.join(output_dir, f'Step02_PLSR_{cv_repeat}x10cv_preds.csv'), index_col=0)

    # Visualization
    plt.close()
    fig = plt.figure(figsize=(24,10))

    jet_9 = cm.get_cmap('jet', 9)(np.arange(9))
    cm.get_cmap('jet', 9)

    plt.subplot(121)
    sns.set_context('poster',font_scale=0.8)
    score = []
    for i in range(cv_repeat):
        score.append(stats.pearsonr(preds_rep.iloc[:,i].values, report['Y'])[0])
    jnk = pd.DataFrame(index = range(cv_repeat), 
                        columns = ['score','iteration'])
    jnk.loc[:,'iteration'] = list(range(cv_repeat))
    jnk.loc[:,'score'] = score
    sns.histplot(score, bins=15, kde=True, color=RdBu_200[185], alpha=0.7)
    sns.rugplot(score, color=RdBu_200[185], lw =3)

    sns.despine()
    plt.xlim(0.35, 0.45)
    plt.xticks([0.35,0.37,0.39,0.41,0.43,0.45], fontsize=20)
    plt.yticks(fontsize=20)
    ylim = plt.ylim()
    plt.vlines(np.median(score), ylim[0], ylim[1], linestyle='-',
            color=jet_9[0], linewidth=4)
    plt.text(np.median(score)+0.001, 16,'Median=%s'%(float('%.2g' % np.median(score))), 
            fontsize=30)
    plt.ylabel('Density', fontsize=25)
    plt.xlabel('R between predicted and actual IDP', fontsize=25)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    median_score = np.median(score)
    median_index = np.argmin(np.abs(np.array(score) - median_score)) 
    print(f"median model index: {median_index}")

    plt.subplot(122)
    sns.set_context('paper', font_scale=0.8)
    sns.set_style('ticks')
    sns.regplot(x=report['Y'], y=preds_rep.iloc[:,median_index], fit_reg=True, 
                scatter_kws={'s': 50, 'linewidths': 1, 'color':RdBu_Hex[185]},
                line_kws={'linewidth':4, 'color': RdBu_Hex[195]}) # 'color':'b'
    sns.despine()
    plt.xlim(0,1.0,)
    plt.ylim(-0.12,1.1)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('Actual IDP', fontsize=25)
    plt.ylabel('Predicted IDP', fontsize=25)
    r = np.float16(stats.pearsonr(preds_rep.iloc[:,median_index], report['Y'])[0])
    plt.text(0.05, 0.92,'$r=%s, p_{sa}<0.01$'%(float('%.2g' % r)), fontsize=30)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)

    # Brainsmash
    f = f'{output_dir}BrainSmash/Step02_BrainSmash2PLSR_model{median_index}_preds.csv'
    csv_cur = pd.read_csv(f, index_col=0)
    Final_perm_cur = pd.Series([stats.pearsonr(report['Y'], csv_cur.iloc[:, x])[0] for x in range(len(csv_cur.T))])

    a = plt.axes([.77, .14, 0.13, 0.2])
    sns.histplot(Final_perm_cur, bins=15, kde=True, color=RdBu_200[185], alpha=0.7)
    sns.rugplot(Final_perm_cur, color=RdBu_200[185], lw =3)
    sns.despine()
    ylim = plt.ylim()
    plt.vlines(median_score, ylim[0], ylim[1], linestyle='-',
            color=RdBu_200[185], linewidth=4, label='Actural r')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylabel("Density", fontsize=12)
    plt.text(0.16, 2.2,'Actural r', fontsize=12)

    fig_name = f'Step02_PLSR_{cv_repeat}x10cv_r2median_{median_index}.png'
    plt.savefig(os.path.join(output_dir, fig_name),  dpi=500, bbox_inches='tight')
    plt.show()

def perform_CV(i, expression, report, clf, cv_strategy, cv_repeat):
    y = report['Y']
    preds = np.zeros((len(y),))
    # r_values = np.zeros((cv_strategy,))
    sel = KFold(n_splits=cv_strategy, shuffle=True, random_state=(123 * (i + 1)))
    print(f'repeat: {i}/{cv_repeat}')
    for fold, (tr_ind, te_ind) in enumerate(sel.split(expression)):
        tr_set = expression.iloc[tr_ind]
        tr_y = report['Y'].iloc[tr_ind]
        te_set = expression.iloc[te_ind]
        te_y = report['Y'].iloc[te_ind]
        model = clf.fit(tr_set, tr_y)
        predictions = model.predict(te_set).flatten()
        preds[te_ind] = predictions
        # r_values[fold] = stats.pearsonr(te_y, predictions)[0]
    return preds


def Step02_Brainsmash(input_file_name, n_permutations, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir

    from brainsmash.mapgen.sampled import Sampled
    from brainsmash.workbench.geo import volume
    # The coordinates 
    template_atlas = nib.load(f'{data_dir}Input_example/cerebellum_17853_nifti.nii')   
    Cere = template_atlas.get_fdata()[template_atlas.get_fdata() != 0]
    cere_coords = np.where(template_atlas.get_fdata() != 0)
    cere_coords = [(cere_coords[0][x],
                    cere_coords[1][x],
                    cere_coords[2][x]) for x in range(len(cere_coords[0]))]
    sample_coordinate = cere_coords
    Cere_atlas = nib.load(f'{data_dir}{input_file_name}')
    sample_brainmap  = Cere_atlas.get_fdata()[template_atlas.get_fdata() != 0]

    # calculate distance matrix
    from brainsmash.workbench.geo import volume
    coord_file = np.array(sample_coordinate)
    output_dir = f'{output_dir}BrainSmash/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filenames_cere = volume(coord_file, output_dir)

    start = time.time()
    ans = np.zeros([n_permutations, sample_brainmap.shape[0]])
    knn = 2500
    kwargs = {'ns': 50,
                'knn': knn,
                'pv': 70}
    gen = Sampled(x=sample_brainmap, D=filenames_cere['D'], index=filenames_cere['index'], resample=False, n_jobs=50, **kwargs)
    surrogate_maps = gen(n = n_permutations)
    pd.DataFrame(surrogate_maps).to_csv(os.path.join(output_dir,f'Step02_BrainSmash_resample_{n_permutations}.csv'))

def Step02_Brainsmash2FG(n_permutations, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir

    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)
    loclist = report[['suit_x','suit_y','suit_z']].values.astype('float')  
    Grey_file = nib.load(os.path.join(data_dir,"Input_example/HCP_S1200_1003_rfMRI_MSMAll_groupPCA_d4500ROW_zcorr.dconn.nii"))
    axis2 = Grey_file.header.get_axis(1)
    cerebellum_mask_ind = np.append(np.where(axis2.__dict__['_name'] == 'CIFTI_STRUCTURE_CEREBELLUM_LEFT'), 
                                    np.where(axis2.__dict__['_name'] == 'CIFTI_STRUCTURE_CEREBELLUM_RIGHT'))
    cope1_file = nib.load(os.path.join(data_dir,"Input_example/cope1.dtseries.nii"))
    output = np.zeros([91282])
    G = pd.read_csv(f'{output_dir}BrainSmash/Step02_BrainSmash_resample_{n_permutations}.csv',index_col=0)
    for i in tqdm(range(len(G))):
        output[cerebellum_mask_ind] = G.loc[i,:]
        img = nib.cifti2.Cifti2Image(output.reshape([1,-1]), nib.cifti2.Cifti2Header(cope1_file.header.matrix))
        dscalar_name = f'{output_dir}BrainSmash/BrainSmash_{i+1}.dscalar.nii'
        nifti_name = f'{output_dir}BrainSmash/BrainSmash_{i+1}_nifti.nii'
        if not os.path.exists(nifti_name):
            img.to_filename(dscalar_name)
            cmd = f"wb_command -cifti-separate {dscalar_name} COLUMN -volume-all {nifti_name}"
            subprocess.check_output(cmd, shell=True)

    nifti_name = f'{output_dir}BrainSmash/BrainSmash_*_nifti.nii'
    xpsheets = sorted(glob(nifti_name))
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    tasks = [[x] for x in sorted(glob(nifti_name))]
    y_gradient_list = pool.starmap(GetFG, tasks)
    y_gradient = pd.concat([pd.DataFrame(y_gradient_list[x]) for x in range(len(y_gradient_list))], axis=1)
    filename = f'{output_dir}BrainSmash/BrainSmash_Rad4mm_sample317_from17853_resample_{y_gradient.shape[1]}.csv'
    pd.DataFrame(y_gradient).to_csv(filename) 

def Step02_Brainsmash2PLSR(n_components, median_index, n_permutations, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir

    expression = pd.read_csv(os.path.join(output_dir, 'Step01_Gene_expression.csv'), index_col=0)
    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)
    report.index = expression.index
    FG_Perm = pd.read_csv(f'{output_dir}BrainSmash/BrainSmash_Rad4mm_sample317_from17853_resample_{n_permutations}.csv', index_col=0)
    FG_Perm.index = report.index
    FG_Perm = FG_Perm.dropna(axis = 1, how = 'any').iloc[:, 0:10000]
    # FG_Perm_cur = FG_Perm.iloc[:, (PTime_each-1000):PTime_each]

    clf = PLSRegression(n_components)
    cv = 10 # 10 folds
    n_jobs = -1  # 使用所有可用的 CPU 核心
    results = Parallel(n_jobs=n_jobs)(
        delayed(Brainsmash_2PLSR)(perm, FG_Perm, expression, report, n_components, cv, clf, median_index)
        for perm in range(len(FG_Perm.T))
    )
    preds = pd.concat([res[0] for res in results], axis=1)
    score = pd.DataFrame([res[1] for res in results]).T

    output_dir = f'{output_dir}BrainSmash/'
    preds_name = f'{output_dir}/Step02_BrainSmash2PLSR_model{median_index}_preds.csv'
    score_name = f'{output_dir}/Step02_BrainSmash2PLSR_model{median_index}_r_1fold.csv'
    preds.to_csv(preds_name)
    score.to_csv(score_name)
    return preds_name, score_name

def Brainsmash_2PLSR(perm, FG_Perm, expression, report, n_components, cv, clf, median_index):
    print(f'repeat: {perm}/{len(FG_Perm.T)}')
    y = FG_Perm.iloc[:, perm]
    preds_cur = pd.DataFrame(np.zeros((len(y), 1)), index=expression.index)
    r_cur_shuffleY = []
    r_cur_originalY = []
    
    sel = model_selection.KFold(n_splits=cv, shuffle=True, random_state=(123 * (median_index+1)))
    for tr_ind, te_ind in sel.split(expression):
        tr_set, tr_y = expression.iloc[tr_ind], y.iloc[tr_ind]
        te_set, te_y_shuffle, te_y_original = expression.iloc[te_ind], y.iloc[te_ind], report['Y'].iloc[te_ind]
        mod = clf.fit(tr_set, tr_y)
        preds_te = mod.predict(te_set)
        preds_cur.iloc[te_ind, 0] = preds_te.ravel()
        r_cur_shuffleY.append(stats.pearsonr(te_y_shuffle, preds_te.ravel())[0])
        r_cur_originalY.append(stats.pearsonr(te_y_original, preds_te.ravel())[0])

    score = [
        stats.pearsonr(preds_cur.iloc[:, 0], y)[0],  # avg2r
        stats.pearsonr(preds_cur.iloc[:, 0], report['Y'])[0],  # r2median_originalY
        stats.pearsonr(preds_cur.iloc[:, 0], y)[0],  # r2median_shuffleY
        np.mean(r_cur_originalY),  # r_1fold_originalY
        np.mean(r_cur_shuffleY)  # r_1fold_shuffleY
    ]

    return preds_cur, score

def disbtw(surf, aa):
    surf= np.array(surf)
    aa = np.expand_dims(aa,axis=0)
    dis = np.sqrt(((surf-aa)**2).sum(axis=1))
    return dis

def Grid2world(GridList, affine):
    """
    trans the grid coords into the corresponding MNI coords in MNI space
    GridList shape: N*3, affine shape 4*4
    """
    GridList = np.array(GridList)
    length = len(GridList)
    affine = np.array(affine)
    Grid_coord = np.concatenate([GridList, np.ones([length, 1])], axis=1)
    MNI_coord = Grid_coord.dot(affine.T)      
    MNI_coord = MNI_coord.round(2)
    return MNI_coord[:, :3]

def GetFG(sheet, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)
    loclist = report[['suit_x','suit_y','suit_z']].values.astype('float')  
    y_gradient = np.zeros([len(loclist), 1])
    fid = int(sheet.split('/')[-1].split('_')[-2]) # get functional gradient id
    # print(fid)
    gradient_file = nib.load(sheet)
    gradient = gradient_file.get_fdata()
    affine = gradient_file.affine
    grid_list = np.nonzero(gradient)   
    grad_value = gradient[grid_list] 
    rrlist = []
    MNI_loc = Grid2world(np.array(grid_list).T, affine)
    for loc in loclist:
        dis = disbtw(MNI_loc, loc)
        tmp = grad_value[dis<4].mean()
        rrlist.append(tmp)
        # print(np.array(rrlist))
    y_gradient = rrlist
    return y_gradient



def Step03_GCIsig(n_components, illustrative, n_permutations, data_dir=None, output_dir=None):
    if data_dir is None:
        data_dir = global_settings.data_dir
    if output_dir is None:
        output_dir = global_settings.output_dir
    expression = pd.read_csv(os.path.join(output_dir, 'Step01_Gene_expression.csv'), index_col=0)
    report = pd.read_csv(os.path.join(output_dir, 'Step01_Sample_info.csv'), index_col=0)
    report.index = expression.index

    final_outputs = {}
    y = report['Y']
    clf = PLSRegression(n_components)
    print('running final model')
    mod = clf.fit(expression, y)
    final_outputs.update({'final_model': mod})
    scr = mod.score(expression, y) 
    print('final model fit r = ',scr)

    f_betas = mod.coef_
    final_outputs.update({'betas': f_betas})
    final_outputs.update({'r': scr})

    if illustrative:
        plt.close()
        sns.regplot(x=mod.predict(expression), y=y)
        plt.xlabel('Model predicted IDP')
        plt.ylabel('Actual IDP')
        plt.show()
        plt.close()

        # fig = plt.figure(figsize=(15,10))
        # sns.set_context('poster',font_scale=0.8)

        # # Major
        # sns.histplot(f_betas, bins=20, kde=False, color=RdBu_200[185], alpha=0.7)
        # sns.despine()
        # plt.xlim(-0.0018,0.0016)
        # plt.xticks([-0.0015,-0.00075,0,0.00075,0.0015], fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.xlabel('PLSR beta', fontsize=25)
        # plt.ylabel('Number', fontsize=25)

        # # left
        # al = plt.axes([.18, .2, 0.15, 0.4])
        # sns.histplot(f_betas, bins=20, kde=False, color=RdBu_200[185], alpha=0.7)
        # plt.xlim(-0.0020,-0.001)
        # plt.xticks([-0.002,-0.0015,-0.001], fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.ylim(0,25)
        # xval = sorted(f_betas)[50]
        # plt.plot([xval, xval], [0, 100], linewidth=4, color=RdBu_200[185])

        # # right
        # ar = plt.axes([.7, .2, 0.15, 0.4])
        # sns.histplot(f_betas, bins=20, kde=False, color=RdBu_200[185], alpha=0.7)
        # plt.xlim(0.0009,0.0016)
        # plt.xticks([0.0010,0.0013,0.0016], fontsize=20)
        # plt.yticks(fontsize=20)
        # plt.ylim(0,25)
        # xval = sorted(f_betas)[-50]
        # plt.plot([xval, xval], [0, 100], linewidth=4, color=RdBu_200[185])


        # plt.savefig(os.path.join(output_dir,'Step03_Finalmodel_GCI.png'),
        #             dpi=700/24*15, bbox_inches='tight')
        # plt.show()
    df = pd.DataFrame(f_betas.T)
    df.index = expression.columns
    df.columns = ['beta']
    gene_table = df.sort_values(ascending=False, by = 'beta')
    # gene_table.to_csv(os.path.join(output_dir,'Step03_GCI_BETA.csv'))


    FG_Perm = pd.read_csv(f'{output_dir}BrainSmash/BrainSmash_Rad4mm_sample317_from17853_resample_{n_permutations}.csv', index_col=0)
    FG_Perm.index = report.index
    FG_Perm = FG_Perm.dropna(axis = 1, how = 'any').iloc[:, 0:10000]
    clf = PLSRegression(n_components)

    n_jobs = -1  # Use all available CPU cores
    results = Parallel(n_jobs=n_jobs)(
        delayed(GCIsig)(perm, FG_Perm, expression, report, n_components,clf)
        for perm in range(len(FG_Perm.T))
    )
    
    score = pd.DataFrame([res[0] for res in results]).T
    BETA = pd.concat([pd.Series(res[1]) for res in results], axis=1)
    score_name = f'{output_dir}BrainSmash/Step03_GCI_score_perm.csv'
    score.to_csv(score_name)
    BETA_name = f'{output_dir}BrainSmash/Step03_GCI_BETA_perm.csv'
    BETA.to_csv(BETA_name)


    Beta_perm= pd.read_csv(BETA_name, index_col=0)
    Beta_perm.index = expression.columns
    Beta_perm.columns = range(len(Beta_perm.T))

    gene_table['p_perm'] = [abs(Beta_perm.loc[x,:]).ge(abs(gene_table.loc[x,'beta'])).sum()/len(Beta_perm.T) for x in gene_table.index]
    # 1/10001 > (0.05/15624), so we just use the value larger than the tru divided by the permutation times
    gene_table['p_perm_fdr0.05'] = gene_table['p_perm']<(0.05/15624)
    print('GCIsig number')
    print(len(gene_table[gene_table['p_perm_fdr0.05']==True]))
    gene_table['type'] = ''
    for x in gene_table.index: 
        if gene_table.loc[x,'beta'] > 0:
            gene_table.loc[x,'type'] = 'pos'
            if (~gene_table.loc[x,'p_perm_fdr0.05']):
                gene_table.loc[x,'type'] = 'nosig'
        if gene_table.loc[x,'beta'] < 0:
            gene_table.loc[x,'type'] = 'neg'
            if (~gene_table.loc[x,'p_perm_fdr0.05']):
                gene_table.loc[x,'type'] = 'nosig'
    gene_table['color'] = ''
    for x in gene_table.index: 
        if gene_table.loc[x,'beta'] > 0:
            gene_table.loc[x,'color'] = Hex_2[1]
            if (~gene_table.loc[x,'p_perm_fdr0.05']):
                gene_table.loc[x,'color'] = "grey"
        if gene_table.loc[x,'beta'] < 0:
            gene_table.loc[x,'color'] = Hex_1[2]
            if (~gene_table.loc[x,'p_perm_fdr0.05']):
                gene_table.loc[x,'color'] = "grey"
    gene_table_file = os.path.join(output_dir,'Step03_GCIsig.csv')
    gene_table.to_csv(gene_table_file)
    # return gene_table_file

def GCIsig(perm, FG_Perm, expression, report, n_components, clf):
    print(f'repeat: {perm}/{len(FG_Perm.T)}')
    y = FG_Perm.iloc[:, perm]
    
    mod = clf.fit(expression, y)
    score = abs(stats.pearsonr(mod.predict(expression).flatten(), y)[0])
    BETA = mod.coef_.flatten()
    return score, BETA

    