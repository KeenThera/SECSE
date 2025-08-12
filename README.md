# **SECSE**

----------------------------

### SECSE: _**S**ystemic **E**volutionary **C**hemical **S**pace **E**xplorer_

![plot](docs/platform.jpg)

Chemical space exploration is a major task of the hit-finding process during the pursuit of novel chemical entities.
Compared with other screening technologies, computational _de novo_ design has become a popular approach to overcome the
limitation of current chemical libraries. Here, we reported a _de novo_ design platform named systemic evolutionary
chemical space explorer (SECSE). The platform was conceptually inspired by fragment-based drug design, that miniaturized
a “lego-building” process within the pocket of a certain target. The key to virtual hits generation was then turned into
a computational search problem. To enhance search and optimization, human intelligence and deep learning were
integrated. SECSE has the potential in finding novel and diverse small molecules that are attractive starting points for
further validation.

### Tutorials and Usage

----------------------------

1. Setting up dependencies  
   python ~=3.9, perl ~=5.32
    ```bash
    conda create --name secse -c conda-forge parallel tqdm biopandas openbabel chemprop xlrd=2 pandarallel rdkit=2024.09.1 loguru tensorboard
    conda activate secse
   ```
2. Installing from source
    ```bash
    git clone https://github.com/KeenThera/SECSE.git 
   ```
3. Setting Environment Variables  
   `export SECSE=/absolute/path/to/SECSE`  
   I'm using AutoDock Vina for docking:
   [(download here)](https://github.com/ccsb-scripps/AutoDock-Vina/releases)  
   `export VINA=/absolute/path/to/AutoDockVINA`  
   I'm using AutoDock GPU: (adgpu-v1.5.3_linux_ocl_128wi)
   [(download here)](https://github.com/ccsb-scripps/AutoDock-GPU/releases)  
   `export AUTODOCK_GPU=/absolute/path/to/AutoDockGPU`  
   I'm using [Gilde](https://www.schrodinger.com/products/glide) for docking (additional installation & license
   required):  
   `export SCHRODINGER=/absolute/path/to/SCHRODINGER`  
   I'm using [Uni-Dock](https://github.com/dptech-corp/Uni-Dock) for docking (need GPU):  
   [compile from Uni-Dock source code](https://github.com/dptech-corp/Uni-Dock/tree/main/unidock#building-from-source) (recommand), or [download here](https://github.com/dptech-corp/Uni-Dock/releases/download/1.1.0/unidock-1.1.0-cuda120-linux-x86_64) and add `export UNIDOCK=/absolute/path/to/UNIDOCK`
4. Giving execution permissions to the SECSE directory  
   `chmod -R +x /absolute/path/to/SECSE`
5. Input fragments: a tab separated _.smi_ file without header. See demo [here](demo/demo_1020.smi).
6. Parameters in config file:

   [general]

    - _project_code_, project identifier, which will be prefixed to each generated molecule ID, type=str
    - _workdir_, working directory, create if not exists, otherwise overwrite, type=str
    - _fragments_, file path to seed fragments, smi format, type=str
    - _num_per_gen_, number of molecules generated each generation, type=int
    - _seed_per_gen_, number of selected seed molecules per generation, default=1000, type=int
    - _start_gen_, number of staring generation, if you want to resume the generation, please specify the 'start_gen' as
      the number corresponding to the last **completed generation** in your previous run, default=0, type=int
    - _num_gen_, number of growing generations, the final generation number will be the sum of start_gen and num_gen,
      type=int
    
    - _cpu_, number of max invoke CPUs, type=int
    - _gpu_, number of max invoke GPU for AutoDock GPU, type=int
    - _rule_db_, path to customized rule in json format, input 0 if use default rule, default=0

   [docking]
    - _docking_program_, name of docking program, AutoDock-Vina (input vina) or AutoDock-GPU (input autodock-gpu) or
      Glide (input glide) , default=vina, type=str
    - _target_, protein PDBQT if use AutoDock Vina; grid map files descriptor fld file if AutoDock GPU; Grid file if
      choose Glide, type=str
    - _RMSD_, docking pose RMSD cutoff between children and parent, default=2, type=float
    - _delta_score_, decreased docking score cutoff between children and parent, default=-1.0, type=float
    - _score_cutoff_, default=-9, type=float

   Parameters when docking by AutoDock Vina:

    - _x_, Docking box x, type=float
    - _y_, Docking box y, type=float
    - _z_, Docking box z, type=float
    - _box_size_x_, Docking box size x, default=20, type=float
    - _box_size_y_, Docking box size y, default=20, type=float
    - _box_size_z_, Docking box size z, default=20, type=float

   [prediction]

    - _mode_, mode of deep learning modeling, 0: not use, 1: modeling per generation, 2: modeling overall after all the
      generation, default=0, type=int
    - _dl_per_gen_, top N predicted molecules for docking, default=100, type=int
    - _dl_score_cutoff_, default=-9, type=float

   [properties]

    - _mw_, molecular weights cutoff, default=450, type=int
    - _logp_lower_, minimum of logP, default=0.5, type=float
    - _logp_upper_, maximum of logP, default=7, type=float
    - _chiral_center_, maximum of chiral center,default=2, type=int
    - _heteroatom_ratio_, maximum of heteroatom ratio, default=0.35, type=float
    - _rdkit_rotatable_bound_num_, maximum of rotatable bound calculated from
      rdkit.rdMolDescriptors.CalcNumRotatableBonds, default=5, type=int
    - _keen_rotatable_bound_num_, maximum of rotatable bound defined by KEEN (
      SMARTS: "[C^3!D1;!$(C(F)(F)F)]-!@[!Br!F!Cl!I!H3&!$(*#*)!D1;!$([!Br!F!Cl!I](F)(F)F)]"), default=3, type=int
    - _rigid_body_num_, maximum of rigid body defined by KEEN (
      SMARTS: "[C^3!D1;!$(C(F)(F)F);!R;!$(C=O(N));!$(NC(=O));!$(C(=O)O);!$(C(=O)O)]-!@[!Br!F!Cl!I!H3&!$(*#*)!
      D1;!$([!Br!F!Cl!I](F)(F)F);!R;!$(C=O([N,O]));!$(NC(=O));!$(C(=O)O)]"), default=2, type=int
    - _hbd_, maximum of hydrogen bond donor calculated by rdkit.rdMolDescriptors.CalcNumHBD, default=5, type=int
    - _hba_, maximum of hydrogen bond acceptor calculated by rdkit.rdMolDescriptors.CalcNumHBA, default=10, type=int
    - _tpsa_, maximum of topological polar surface area calculated by rdkit.Chem.Descriptors.TPSA, default=200,
      type=float
    - _lipinski_violation_, maximum of violation of Lipinski rule of five calculated by RDKit, default=1, default=1,
      type=int
    - _qed_, QED (calculated by rdkit.Chem.QED.qed) cutoff value, default=0.5, type=float
    - _max_ring_size_, maximum of ring size, default=7, type=int
    - _max_ring_system_size_, maximum of ring system member size in one ring system, default=3, type=int
    - _ring_system_count_, maximum of seperated ring system count, default=4, type=int
    - _bridged_site_count_, maximum of bridged ring site count, default=2, type=int
    - _spiro_site_count_, maximum of spiro ring site count, default=1, type=int
    - _fused_site_count_, maximum of fused ring site count, default=3, type=int
    - _rdkit_sa_score_, synthetic accessibility score (calculated by RDKit) cutoff, default=5, type=float
    - _substructure_filter_, files containing the customized unwanted substructure SMARTS in "*.xls" format, set the
      value to 0 if you do not have any additional unwanted substructure. PANIS already includes as default. The file
      should include columns for **`Pattern`**,  **`ID`**, and **`Max`**, where the **`ID`** should be unique for each SMARTS. You can
      refer to the example file [subtructure_filter_demo.xls](demo/subtructure_filter_demo.xls), default=0, type=string

   Config file of a demo case [phgdh_demo_vina.ini](demo/phgdh_demo_vina.ini)  
   Customized rule json template [rules.json](demo/rules.json). Rule ID should be in the form G-001-XXXX, like
   G-001-0001, G-001-0002, G-001-0003 ...

7. Run SECSE  
   `python $SECSE/run_secse.py --config /absolute/path/to/config`  
   Please input the **absolute path** of the config file here.
8. Output files
    - merged_docked_best_timestamp_with_grow_path.csv: selected molecules and growing path
    - selected.sdf: 3D conformers of all selected molecules

### Dependencies

-------
GNU Parallel installation

- CentOS / RHEL  
  `sudo yum install parallel`
- Ubuntu / Debian  
  `sudo apt-get install parallel`
- From source: https://www.gnu.org/software/parallel/

python ~=3.12, perl ~=5.32

numpy~=1.26.4, pandas~=2.2.2, xlrd～=2.0.1, pandarallel~=1.6.5, tqdm~=4.67.0, biopandas~=0.5.1, openbabel~=3.1.1, rdkit~
=2024.09, chemprop~=2.1, pytorch~=2.5.1+cu117, tensorboard~=2.18.0

Linux server with CPUs only also works.

### Citation

-------
Lu, C.; Liu, S.; Shi, W.; Yu, J.; Zhou, Z.; Zhang, X.; Lu, X.; Cai, F.; Xia, N.; Wang, Y. Systemic Evolutionary Chemical
Space Exploration For Drug Discovery. J Cheminform 14, 19 (2022).

https://doi.org/10.1186/s13321-022-00598-4

### License

-------
SECSE is released under [Apache License, Version 2.0](LICENSE.txt).

The project is being actively developed, if you have any questions or suggestions, please contact:
wang_yikai@keenthera.com or luchong121@outlook.com
