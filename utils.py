import pandas as pd
from OT_deep_score_src.evaluate_utilities import measure_clevage_acc
from OT_deep_score_src.general_utilities import Model_task
from train_and_predict_scripts.utilities import ensemble_predict

def ensemble_predict_wrapper():
    ensemble_predict(
            ensemble_components_file_path_and_name_list=[
               "/home/bcrlab/barsapi1/CRISPR-Bulge-V/files/bulges/1_folds/Transfer_Learning/ensemble/epochs_10/5_revision_exclude_NewGUIDEseq_ensemble_{}_continue_from_change_seq/read_ts_0/cleavage_models/aligned/FullGUIDEseq/regression/c_3/ln_x_plus_one_trans/model_fold_0".format(i) for i in range(5)],
            dataset_df="files/datasets/NewGUIDEseq/include_on_targets/NewGUIDEseq_CR_Lazzarotto_2020_dataset.csv")
    
def measure_clevage_acc_wrapper(prediction_file_name):
    return measure_clevage_acc(
        prediction_file_name=prediction_file_name,
        models=["pred_averege_ensemble", ],
        only_bulges=True,
        only_mismatches=False,
        model_task=Model_task.REGRESSION_TASK
    )