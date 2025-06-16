#### This is the readme file that explain how to use the methods.

**Note**: To run these scripts, you would have to modify the parameters in the scripts. All `json` or `txt` files just indicate a format, please use customized files when running the scripts

In all experiments, `Mock1Abundance.csv`, `Mock2Abundance.tsv` and `Mock3Abundance.csv` indicates the abundance of **3** Mocks in the experimemts. 

The `/seq/folder` would contain the sequences files with the species name in the file name.

## Run classification

### aodp

The file `run_aodp.sh` is used for running the mapping pipeline. The file `analyze_result_part_all.sh` is used for analyze the results and generating the output json file.

### minimap2
The `map_pipeline.sh` is used for mapping the sequences to the reference. The `check_result_part.py` is used to read the file output by minimap2 and generate the file into json format. Which is used for later analysis.

### kraken2
The `map_pipeline.sh` is running kraken2 against the customized database. In the result, some sequences may be mapped to a level that higher than species, so there are two scripts `find_high_level_taxonomy.py` which goes over the output files and find the ncbi taxo id that are taxonomy ranks that higher than species, `write_taxonomy_high_level.py` which write the taxonomy rank of those ids to `ranks_app.json`. Then running the `analyze_result_part_all.sh` would be able to get the json format result.

### EMU

The `build_taxonomy_tsv.py` is used to generate the file `taxonomy.tsv`, which is later used for building the customized database.

```
emu build-database <db_name> --sequences <database.fasta> --seq2tax <seq2taxid.map> --taxonomy-list <taxonomy.tsv>
```
Then run `run_emu.sh` to do the mapping. And run `analyze_result_part_all.sh` to genearte the json format result file.

### DNABarcoder

The `write_Data.py` is used to generate the `ref.classification` file, which is later used for running DNABarcoder.
Then run the `run_dnabarcoder_part.sh` to run DNABarcoder. In this method, the input is conbined fasta file, which means that all the sequences from different species are in one single fasta file, when we do this experiment, we combined the sequences with species in alphabetical order, which is convenient for evaluation in the future.

Then run `check_result_all.sh` and get the json output files


### how to select the species that used for testing overfitting on vuthuyduong_cnn and MycoAI methods

You need to build a Docker image before running the two scripts `find_species.sh` and `summarize_var.sh`. 

Run `find_species.sh` using the command below.
```
./find_species.sh /fasta/files/folder /output/folder
```
Run `summarize_var.sh` using the command below. The parameter is the `/output/folder` in the last command.
```
./summarize_var.sh /output/folder
```

### vuthuyduong_cnn

In this method, the train set or test set is conbined fasta file, which means that all the sequences from different species are in one single fasta file, when we do this experiment, we combined the sequences with species in alphabetical order, which is convenient for evaluation in the future.

You will need to run `write_data.py` to generate the file `train_taxa.txt`, which stores the taxonomy ranks for each sequences in the train set, this script would also combine train sequences from different species to one fasta file.

Run the python command to train the model, you may need to modify the parameters in `trainCNN.py` before that.
```
python3 trainCNN.py -i train_sequences.fasta -c train_taxa.txt -p 6
```

File `classify_part.sh` is used for classify the sequences in the test set.

File `check_result_part.sh` is used for analyzing the result and generating the json files.

## test overfitting

Run `classify_test_overfitting.sh` to classify sequences with model generated during the training process after certain numbers of epochs.

Run `check_result_testing_overfit.sh` to analyze the classification result of specific species on the model generated during the training process. you can modify line 11 together with file `Mock1Abundance_test_overfit.csv` to specific on other species.

Ths file `plot_json_trends.py` is used for plot the trends of corrected classified sequences number for selected species during training. Can modify the command below when running.

```
python3 plot_json_trends.py /json/file/folder vuthuyduong Qmin17 out_folder

```


### MycoAI

File `write_data.py` is used for read the train set and write all of them to a fasta with each sequence has id and description format meet the requirement of MycoAI input.

If you want to save the model during the training process, use `train_save_middle.sh`, or use `train.sh` if you only want the final model. you can modify the script to add `--base_arch_type CNN` or not to choose using `BERT` model or `CNN` model.

You can modify line 53 to 59 in `train_model_with_mid.py` to determine how many models you would like to sane during the training process.

`classify_part.sh` is used for classify the sequences using trained model. `evaluate_json_part_all.sh` is used for analyzing the results and generating the json file.

## test overfitting

`classify_test_overfit.sh` is used for classify the sequences of the species used for testing overfitting.
`evaluate_json_test_overfit.sh` is used for analyzing the results classified by the model saved during the training process and generating the json file.

Ths file `plot_json_trends.py` is used for plot the trends of corrected classified sequences number for selected species during training. Can modify the command below when running.

```
python3 plot_json_trends.py /json/files/folder cnn Qmin17 /output/folder

```


## plot the figures

See the example command in `draw_recall.sh`, `draw_precision.sh` and `draw_f1.sh`. The `gsd_path.txt` shows an example, in each line, there is a folder path that contains json file in it. The order of the folders would be the order of the plots of different methods.

### JSON file format
In this paper, we used `json` format to show the results, for each json file, it contains a dictionary generated by python script, the key is a species used for analysis, the value is a sub-dictionary. In the sub-dictionary, the key is the taxonomy rank of a species (in the result of kraken2, it may be a taxonomy rank in a higher level), the value is how many sequences classified into this species (or higher taxonomy rank in the result of kraken2).



