###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
G = 6.674e-8 # 万有引力常数
###################################################
import math
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import matplotlib.pyplot as plt
import h5py
import numpy as np

font_path = '/home/cjq/Downloads/hybridCONe_plot/hybridCONe_plot/TimesNewRoman.ttf'
font_prop_label = fm.FontProperties(fname=font_path,size=18)

def categorize_particles_by_y(snapshot_file, particle_type):
    """
    根据x坐标的正负将粒子分为两组，并提取它们的物理量
    
    参数:
    snapshot_file: Gadget4输出的HDF5文件路径
    particle_type: 要分析的粒子类型，默认为PartType4(恒星)
    
    返回:
    一个字典，包含两组粒子的所有物理量
    """
    # 初始化存储结果的字典
    results = {
        'star1': {},  # x < 0 的粒子
        'star2': {}   # x >= 0 的粒子
    }
    
    with h5py.File(snapshot_file, 'r') as f:
        # 1. 读取坐标并创建分类掩码
        coords = f[particle_type]['Coordinates'][:]
        x_coords = coords[:, 0]
        
        # 创建分类掩码
        mask_star1 = x_coords < 0.5*R_to_solar   # x < 0
        mask_star2 = x_coords >= 0.5*R_to_solar  # x >= 0
        

        
        # 2. 提取并存储坐标
        results['star1']['Coordinates'] = coords[mask_star1]
        results['star2']['Coordinates'] = coords[mask_star2]
        
        # 3. 提取并存储速度
        if 'Velocities' in f[particle_type]:
            velocities = f[particle_type]['Velocities'][:]
            results['star1']['Velocities'] = velocities[mask_star1]
            results['star2']['Velocities'] = velocities[mask_star2]
        
        # 4. 提取并存储内能 (对于气体粒子PartType0)
        if particle_type == 'PartType0' and 'InternalEnergy' in f[particle_type]:
            internal_energy = f[particle_type]['InternalEnergy'][:]
            results['star1']['InternalEnergy'] = internal_energy[mask_star1]
            results['star2']['InternalEnergy'] = internal_energy[mask_star2]
        
        # 5. 提取并存储质量
        if 'Masses' in f[particle_type]:
            masses = f[particle_type]['Masses'][:]
            results['star1']['Masses'] = masses[mask_star1]
            results['star2']['Masses'] = masses[mask_star2]
        else:
            # 如果粒子是等质量的，从Header中获取质量
            mass_table = f['Header'].attrs['MassTable']
            type_index = int(particle_type[-1])  # 提取粒子类型编号
            mass_per_particle = mass_table[type_index]
            
            # 为每组粒子创建质量数组
            results['star1']['Masses'] = np.full(np.sum(mask_star1), mass_per_particle)
            results['star2']['Masses'] = np.full(np.sum(mask_star2), mass_per_particle)
        
        # 6. 提取并存储粒子ID (用于跟踪)
        if 'ParticleIDs' in f[particle_type]:
            particle_ids = f[particle_type]['ParticleIDs'][:]
            results['star1']['ParticleIDs'] = particle_ids[mask_star1]
            results['star2']['ParticleIDs'] = particle_ids[mask_star2]
            
            
        if 'Potential' in f[particle_type]:
            particle_ids = f[particle_type]['Potential'][:]
            results['star1']['Potential'] = particle_ids[mask_star1]
            results['star2']['Potential'] = particle_ids[mask_star2]
            
    return results

def tarck_particles_by_ids(snapshot_file,target_ids_star1,target_ids_star2,particle_type):
    results = {
        'star1': {},  # x < 0 的粒子
        'star2': {}   # x >= 0 的粒子
    }
    with h5py.File(snapshot_file, 'r') as f:
        all_ids = f[particle_type]['ParticleIDs'][:]
        mask1 = np.isin(all_ids,target_ids_star1)
        mask2 = np.isin(all_ids,target_ids_star2)
        
        found_ids_star1 = all_ids[mask1]
        found_ids_star2 = all_ids[mask2]
        
        if len(found_ids_star1) < len(target_ids_star1):
            missing_ids = set(target_ids_star1) - set(found_ids_star1)
            print(f"Warning:Missing{len(missing_ids)} particles in donor")
        
        if len(found_ids_star2) < len(target_ids_star2):
            missing_ids = set(target_ids_star2) - set(found_ids_star2)
            print(f"Warning:Missing{len(missing_ids)} particles in accretor")
        
        coords = f[particle_type]['Coordinates'][:]
        velocities = f[particle_type]['Velocities'][:]
        internal_energy = f[particle_type]['InternalEnergy'][:]
        masses = f[particle_type]['Masses'][:]
        particle_ids = f[particle_type]['Potential'][:]
        
        results['star1']['Coordinates'] = coords[mask1]
        results['star2']['Coordinates'] = coords[mask2]
        results['star1']['Velocities'] = velocities[mask1]
        results['star2']['Velocities'] = velocities[mask2]
        results['star1']['InternalEnergy'] = internal_energy[mask1]
        results['star2']['InternalEnergy'] = internal_energy[mask2]
        results['star1']['Masses'] = masses[mask1]
        results['star2']['Masses'] = masses[mask2]
        results['star1']['ParticleIDs'] = particle_ids[mask1]
        results['star2']['ParticleIDs'] = particle_ids[mask2]
        results['star1']['Potential'] = particle_ids[mask1]
        results['star2']['Potential'] = particle_ids[mask2]
            
    return results
            
            
            