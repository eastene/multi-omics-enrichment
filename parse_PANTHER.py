import json

import pandas as pd

# value between 0 (lowest level) and 7 (most abstract level)
anno_level = 9

json_file = "GO_CEL_COMP.json"
json_contents = None

with open(json_file, 'r') as fj:
    json_contents = json.load(fj)

json_top = json_contents["overrepresentation"]
json_out = json_top["group"]


row_dict = {"ID": None, "Level": None, "Pathway": None, "Pathway Size": None, 
            "Expected": None, "Observed": None, "Fold Enrichment": None,
            "p-value (Fisher's)": None}
out_df = pd.DataFrame(columns=row_dict.keys())

for results in json_out:
    res_iter = results["result"]
    if type(res_iter) == dict:
        result = res_iter
        anno_term = result["term"]
        if anno_term["label"] == "UNCLASSIFIED":
            row_dict["Pathway"] = anno_term["label"]
            row_dict["ID"] = "N/A"
            row_dict["Level"] = "N/A"

            row_dict["Pathway Size"] = result["number_in_reference"]

            result_contents = result["input_list"]
            row_dict["Observed"] = result_contents["number_in_list"]
            row_dict["Expected"] = result_contents["expected"]
            row_dict["p-value (Fisher's)"] = result_contents["pValue"]
            row_dict["Fold Enrichment"] = result_contents["fold_enrichment"]
            out_df = out_df.append(row_dict, ignore_index=True)

        elif anno_term["level"] == anno_level:
            row_dict["Level"] = anno_term["level"]
            row_dict["Pathway"] = anno_term["label"]
            row_dict["ID"] = anno_term["id"]

            row_dict["Pathway Size"] = result["number_in_reference"]

            result_contents = result["input_list"]
            row_dict["Observed"] = result_contents["number_in_list"]
            row_dict["Expected"] = result_contents["expected"]
            row_dict["p-value (Fisher's)"] = result_contents["pValue"]
            row_dict["Fold Enrichment"] = result_contents["fold_enrichment"]
            out_df = out_df.append(row_dict, ignore_index=True)

    else:
        for result in res_iter:
            anno_term = result["term"]
            if anno_term["level"] == anno_level:
                row_dict["Level"] = anno_term["level"]
                row_dict["Pathway"] = anno_term["label"]
                row_dict["ID"] = anno_term["id"]

                row_dict["Pathway Size"] = result["number_in_reference"]

                result_contents = result["input_list"]
                row_dict["Observed"] = result_contents["number_in_list"]
                row_dict["Expected"] = result_contents["expected"]
                row_dict["p-value (Fisher's)"] = result_contents["pValue"]
                row_dict["Fold Enrichment"] = result_contents["fold_enrichment"]

                out_df = out_df.append(row_dict, ignore_index=True)

out_df.sort_values(by="p-value (Fisher's)")
