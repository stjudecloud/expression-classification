#!/usr/bin/env python3

import sys
from typing import List
import argparse

import json
import dxpy

def bulk_describe_with_id(file_links: List[str], proj_id: str) -> List:
    """Make api call to get file description in bulk
    API call to /system/findDataObjects with a list of file ids.
    Limited to a project and 1000 results per call
    The file information of interest are size, name, id and folder.
    Args:
        file_links: list of file ids (strings)
        proj_id: project id where described items is found
    Returns:
        list of dictionary containing file infomation
    """

    payload = {
        "id": file_links,
        "describe": {
            "defaultFields": False,
            "fields": {
                "folder": True, "id": True, "name": True, "properties": True,
            }
        },
        "scope": {"project": proj_id},
        "visibility": "either",
        "state": "closed",
    }

    response = dxpy.api.system_find_data_objects(payload)
    results = [x["describe"] for x in response["results"]]

    return results

def chunk_ids(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def main():
    parser = argparse.ArgumentParser()

    # arguments
    parser.add_argument("-p", "--project", help="DNAnexus project identifier")
    parser.add_argument("--ids", help="DNAnexus file ID(s) to search (string)", nargs="+", type=str)
    
    args = parser.parse_args()
    #print(len(args.ids))
    #print(args.ids)
    subsets = chunk_ids(args.ids, 1000)
    results = []
    for s in subsets:
        results = results + bulk_describe_with_id(proj_id=args.project, file_links=s)
    print(json.dumps(results))
    #print(results)
    #print(len(results))

if __name__ == "__main__":
    main()
