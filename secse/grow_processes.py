#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: grow_processes.py
@time: 2021/11/17/13:49
"""
import csv
import shutil
import os
import pandas as pd
import rdkit
from pathlib import Path
from itertools import groupby, islice
from operator import itemgetter
from tempfile import TemporaryDirectory
from evaluate.glide_docking import dock_by_glide
from growing.mutation.mutation import mutation_df
from scoring.ranking import Ranking
from scoring.diversity_score import clustering
from scoring.docking_score_prediction import prepare_files
from scoring.sampling import sample_by_rule_weight
from evaluate.docking import dock_by_py_vina, dock_by_py_autodock_gpu, dock_by_unidock
from evaluate.sdf_output import validate_scored_sdf
from report.grow_path import write_growth
from utilities.load_rules import json_to_DB
from utilities.config import SOURCE_DIR, docking_backend, read_config
from utilities.function_helper import run_command
import time
from loguru import logger

rdkit.RDLogger.DisableLog("rdApp.*")


class Grow(object):
    def __init__(self, generation, mols_smi, workdir, num_per_gen, docking_program, receptor, start_gen, dl_mode,
                 config_path, cpu_num=0, gpu_num=1, rule_db=0, project_code="GEN", x=0, y=0, z=0, box_size_x=0,
                 box_size_y=0, box_size_z=0):

        self.mols_smi = mols_smi
        self.total_generation = int(generation)
        self.workdir = workdir
        self.num_per_gen = num_per_gen
        self.cpu_num = cpu_num
        self.gpu_num = gpu_num

        self.target = receptor
        self.x = x
        self.y = y
        self.z = z
        self.box_size_x = box_size_x
        self.box_size_y = box_size_y
        self.box_size_z = box_size_z

        self.start_gen = start_gen  # record start
        self.gen = start_gen  # generation num for now
        self.docking_program = docking_backend(docking_program)
        self.dl_mode = dl_mode

        self.config_path = config_path

        rule_db = str(rule_db)
        self.rule_json = None
        if rule_db == "0":
            self.rule_db = None
        elif rule_db.endswith("json"):
            self.rule_json = rule_db
            self.rule_db = os.path.join(self.workdir, "rules.db")
        elif rule_db.endswith("db"):
            if not Path(rule_db).is_file():
                raise FileNotFoundError(f"Mutation rule database: {rule_db}")
            self.rule_db = rule_db
        else:
            raise Exception("Please check your input rule file.")
        self.project_code = project_code

        self.lig_sdf = None
        self.winner_df = None
        self.winner_path = None
        self._generation_dir = None
        self._filter_df = None
        self._dock_df = None
        self._sampled_df = None
        self.workdir_now = None

        self._GROWING_STATE_LIST = ["GROWING", "BROKEN", "STOP"]
        self.growing_flag = self._GROWING_STATE_LIST[0]

    def docking_sh(self, step):
        start = time.time()
        os.makedirs(self.workdir_now, exist_ok=True)

        backend = {
            "vina": self.docking_vina,
            "glide": self.docking_glide,
            "autodock-gpu": self.docking_autodock_gpu,
            "unidock": self.docking_unidock,
        }[self.docking_program]
        backend(step)

        # ranking and find top fragments
        self.lig_sdf = os.path.join(self.workdir_now, "docking_outputs_with_score.sdf")
        end = time.time()
        logger.info("Docking time cost: {} min.".format(round((end - start) / 60, 2)))

    def docking_autodock_gpu(self, step):
        logger.info("Step {}: Docking with AutoDock GPU ...".format(step))
        dock_by_py_autodock_gpu(self.workdir_now, self.mols_smi, self.target, self.cpu_num, self.gpu_num)

    def docking_vina(self, step):
        logger.info("Step {}: Docking with Autodock Vina ...".format(step))
        dock_by_py_vina(self.workdir_now, self.mols_smi, self.target, self.cpu_num, self.x, self.y, self.z,
                        self.box_size_x, self.box_size_y, self.box_size_z)

    def docking_glide(self, step):
        logger.info("Step {}: Docking with Glide ...".format(step))
        # set different docking precision for different generation
        if self.gen < 1:
            dock_mode = "SP"
        else:
            dock_mode = "HTVS"
        dock_by_glide(self.workdir_now, self.mols_smi, self.target, self.gen, dock_mode, self.cpu_num)

    def docking_unidock(self, step):
        logger.info("Step {}: Docking with UniDock ...".format(step))
        dock_by_unidock(self.workdir_now, self.mols_smi, self.target, self.cpu_num, self.x, self.y, self.z,
                        self.box_size_x, self.box_size_y, self.box_size_z)

    def ranking_docked_mols(self, step=2):
        logger.info("Step {}: Ranking docked molecules...".format(str(step)))
        ranking = Ranking(sdf=self.lig_sdf, gen=self.gen, config_file=self.config_path)
        if ranking.ranking_flag:
            ranking.docked_df.to_csv(
                os.path.join(self.workdir, "generation_" + str(self.gen), "docked_gen_" + str(self.gen) + ".csv"),
                index=False)
            ranking.tournament_selection()
            # merge mols whose evaluate score below the cutoff
            ranking.mols_score_below_cutoff()
            self.winner_df = ranking.final_df
            # generate smi file
            self.winner_path = os.path.join(self.workdir, "generation_" + str(self.gen),
                                            "best_fragment_gen_" + str(self.gen) + ".smi")
            self.winner_df["id_gen_" + str(self.gen)] = self.winner_df["id_gen_" + str(self.gen)].apply(
                lambda x: x.split("\t")[0])
            self.winner_df[["smiles_gen_" + str(self.gen), "id_gen_" + str(self.gen)]].to_csv(self.winner_path,
                                                                                              sep="\t",
                                                                                              index=False,
                                                                                              quoting=csv.QUOTE_NONE)
        else:
            self.growing_flag = self._GROWING_STATE_LIST[1]
        self.check_growing()

    def dl_pre(self, step):
        logger.info("Step {}.1: Building deep learning models...".format(str(step)))

        train, pre = prepare_files(self.gen, self.workdir, self.dl_mode)
        if pre is None:
            logger.info("Skipping docking score prediction as all molecules have been docked.")
            if self.dl_mode == 2:
                self.dl_mode = 0
            return
        with open(train, newline="") as source:
            # Chemprop needs enough examples to keep at least two in every split.
            enough_training_data = len(list(islice(csv.reader(source), 9))) == 9
        if not enough_training_data:
            logger.info("Skipping docking score prediction: fewer than eight scored training molecules.")
            if self.dl_mode == 2:
                self.dl_mode = 0
            return
        dl_shell = SOURCE_DIR / "scoring" / "chemprop_pre.sh"
        config = read_config(self.config_path)
        dl_select_num = config.get("prediction", "dl_per_gen")
        run_command([dl_shell, self.workdir, train, pre, self.gen, dl_select_num, 22], capture_mode=0)
        # docking top predicted compounds
        self.workdir_now = os.path.join(self.workdir, "generation_{}_pre".format(self.gen))
        self.mols_smi = os.path.join(self.workdir_now, "mols_for_docking_pred.smi")
        self.docking_sh(str(step) + ".2")

        # merge results to the current generation if prediction per generation
        if self.dl_mode == 1:
            self.lig_sdf = os.path.join(self.workdir, "generation_{}".format(self.gen),
                                        "docking_outputs_with_score.sdf")
            docked_sdf = Path(self.lig_sdf)
            prediction_sdf = Path(self.workdir_now) / docked_sdf.name
            with TemporaryDirectory(prefix=".secse-prediction-", dir=docked_sdf.parent) as directory:
                staged = Path(directory) / docked_sdf.name
                with staged.open("wb") as destination:
                    for path in (docked_sdf, prediction_sdf):
                        with path.open("rb") as source:
                            shutil.copyfileobj(source, destination)
                staged.replace(docked_sdf)
            self.workdir_now = os.path.join(self.workdir, "generation_{}".format(self.gen))

    def check_growing(self):
        if self.growing_flag == self._GROWING_STATE_LIST[1]:
            # broken and report generated molecules
            if self.dl_mode == 2:
                self.dl_mode = 0
            write_growth(self.config_path, self.gen - 1, self.dl_mode)
            logger.info("Stopped after generation {} because no candidates remained.", self.gen - 1)
        elif self.growing_flag == self._GROWING_STATE_LIST[2]:
            # regular finsh and stop the program
            write_growth(self.config_path, self.gen, self.dl_mode)
            logger.info("Finished generations {} through {}.", self.start_gen, self.gen)

    def restore_generation(self):
        generation = Path(self.workdir) / f"generation_{self.gen}"
        scores_path = generation / f"docked_gen_{self.gen}.csv"
        seeds_path = generation / f"best_fragment_gen_{self.gen}.smi"
        sdf_path = generation / "docking_outputs_with_score.sdf"
        filter_path = generation / "filter.csv"
        for path in (scores_path, seeds_path, sdf_path, filter_path):
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Cannot resume generation {self.gen}: missing or empty {path}")
        validate_scored_sdf(sdf_path)
        scores = pd.read_csv(scores_path, dtype={"id": str}, keep_default_na=False)
        id_column = f"id_gen_{self.gen}"
        smiles_column = f"smiles_gen_{self.gen}"
        seeds = pd.read_csv(seeds_path, sep="\t", dtype=str, keep_default_na=False)
        if (seeds.empty or not {id_column, smiles_column}.issubset(seeds.columns)
                or not {"id", "smiles", "docking score", "le_ln"}.issubset(scores.columns)):
            raise ValueError(f"Cannot resume generation {self.gen}: invalid seed or ranking columns in {generation}")
        for column in ("docking score", "le_ln"):
            scores[column] = pd.to_numeric(scores[column], errors="raise")
        scores["id"] = scores["id"].str.split("\t", n=1).str[0]
        scores = scores.add_suffix(f"_gen_{self.gen}")
        keys = [smiles_column, id_column]
        ranked_seeds = scores.drop_duplicates(subset=keys).set_index(keys)
        try:
            self.winner_df = ranked_seeds.loc[pd.MultiIndex.from_frame(seeds[keys])].reset_index()
        except KeyError as error:
            raise ValueError(f"Cannot resume: seeds in {seeds_path} are absent from {scores_path}") from error
        self.winner_df = self.winner_df[scores.columns]
        self.winner_path = str(seeds_path)
        self.lig_sdf = str(sdf_path)
        self.workdir_now = str(generation)

    def grow(self):
        logger.info(f"Input fragment file: {self.mols_smi}")
        logger.info(f"Target grid file: {self.target}")
        logger.info(f"Workdir: {self.workdir}")
        logger.info(f"Generation {self.gen} ...")
        step = 2
        if self.start_gen:
            self.restore_generation()
        if self.rule_json is not None and self.total_generation:
            os.makedirs(self.workdir, exist_ok=True)
            json_to_DB(self.rule_json, self.rule_db)
        if not self.start_gen:
            self.workdir_now = os.path.join(self.workdir, "generation_0")
            self.docking_sh(1)
            self.ranking_docked_mols(step)

        # next generations: 1.copy the best mols from last generation as seed; 2.mutation; 3.filter; 4. sampling;
        #                   5.clustering; 6.evaluate; 7.ranking
        for _ in range(self.total_generation):
            self.gen += 1
            logger.info(f"Generation {self.gen} ...")
            self.workdir_now = os.path.join(self.workdir, "generation_" + str(self.gen))
            if os.path.exists(self.workdir_now):
                shutil.rmtree(self.workdir_now)
            os.makedirs(self.workdir_now, exist_ok=True)
            self.winner_df.to_csv(os.path.join(self.workdir_now, "seed_fragments.smi"), sep="\t", index=False,
                                  quoting=csv.QUOTE_NONE)
            # mutation
            logger.info("Step 1: Mutation")

            self._generation_dir = os.path.join(self.workdir_now, "generation_split_by_seed")
            self.winner_df = self.winner_df.reset_index(drop=True)
            header = mutation_df(self.winner_df, self.workdir, self.cpu_num, self.gen, self.rule_db, self.project_code)
            # mutation.csv already contains the first row for each product SMILES.
            generation_path = Path(self.workdir_now)
            shutil.copyfile(generation_path / "mutation.csv", generation_path / "generation.csv")
            split_dir = Path(self._generation_dir)
            split_dir.mkdir()
            with (generation_path / "generation.csv").open(newline="") as source:
                for seed_id, rows in groupby(csv.reader(source), key=itemgetter(1)):
                    if Path(seed_id).name != seed_id or not seed_id:
                        raise ValueError(f"Invalid seed ID for file name: {seed_id!r}")
                    with (split_dir / f"{seed_id}.csv").open("a", newline="") as destination:
                        csv.writer(destination, lineterminator="\n").writerows(rows)
            # filter
            logger.info("Step 2: Applying filter to all mutated molecules.")
            time1 = time.time()
            cmd_filter = [SOURCE_DIR / "growing" / "filter_parallel.sh", self.workdir_now,
                          self.gen, self.config_path, self.cpu_num]
            run_command(cmd_filter)
            time2 = time.time()
            logger.info("Filter runtime: {:.2f} min.".format((time2 - time1) / 60))

            # do not sample or clustering if generated molecules less than wanted size
            try:
                id_columns = {i: str for i, name in enumerate(header) if name.startswith("id_gen_")}
                self._filter_df = pd.read_csv(os.path.join(self.workdir_now, "filter.csv"), header=None,
                                             dtype=id_columns)
            except pd.errors.EmptyDataError:
                self.growing_flag = self._GROWING_STATE_LIST[1]
                logger.info("No molecules met the filter criteria. Please adjust your configuration.")
                self.check_growing()
                return

            self._filter_df.columns = header + ["flag"]
            self._filter_df["type"] = self._filter_df["reaction_id_gen_" + str(self.gen)].apply(
                lambda x: "-".join(x.split("-")[:2]))
            self._filter_df.to_csv(os.path.join(self.workdir_now, "filter.csv"), index=False)
            if self._filter_df.shape[0] <= self.num_per_gen:
                self._dock_df = self._filter_df
                self._dock_df.to_csv(os.path.join(self.workdir_now, "sampled.csv"), index=False)
            else:
                # sampling
                logger.info("Step 3: Sampling")
                self._sampled_df = sample_by_rule_weight(self.gen, self._filter_df, self.workdir_now)
                # self._sampled_df = sample_by_similarity(self.gen, self._filter_df, self.workdir_now, self.num_per_gen)
                logger.info("Step 4: Clustering")
                # clustering
                num_clusters = int(self.num_per_gen / 5) + 1
                self._sampled_df = clustering(self._sampled_df, "smiles_gen_" + str(self.gen), self.gen, self.cpu_num,
                                              num_clusters)

                # sample enough mol
                self._dock_df = self._sampled_df.sort_values("cluster_center_dis_gen_" + str(self.gen)).groupby(
                    "cluster_center_gen_" + str(self.gen)).head(int(self.num_per_gen / num_clusters) + 1)

            # write file for evaluate
            self.mols_smi = os.path.join(self.workdir_now, "mols_for_docking.smi")
            self._dock_df[["smiles_gen_" + str(self.gen), "id_gen_" + str(self.gen)]].to_csv(self.mols_smi, index=False,
                                                                                             header=False, sep="\t")

            # evaluate
            step = 5
            self.docking_sh(step)
            # run deep learning model, when ( dl_mode is 1) & (not all generated compounds were docked)
            if (self.dl_mode == 1) and (self._filter_df.shape[0] > self._dock_df.shape[0]):
                step += 1
                self.dl_pre(step)
            # ranking
            step += 1
            self.ranking_docked_mols(step)
            if self.growing_flag != self._GROWING_STATE_LIST[0]:
                return

        if self.dl_mode == 2 and self.gen > 0:
            step += 1
            self.dl_pre(step)
        elif self.gen == 0:
            self.dl_mode = 0

        self.growing_flag = self._GROWING_STATE_LIST[2]
        self.check_growing()
