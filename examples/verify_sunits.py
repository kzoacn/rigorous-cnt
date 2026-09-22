"""Recheck an exported S-unit record, ignoring its saved success flags."""

import argparse
import json
from pathlib import Path
from rigorous_cnt import verify_sunit_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    args=parser.parse_args()
    result=verify_sunit_record(json.loads(args.record.read_text(encoding='utf-8')))
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if not result['verified']:
        raise SystemExit(1)


if __name__=='__main__':main()
