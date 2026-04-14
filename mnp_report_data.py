#!/usr/bin/env python3
"""Build MNP report summary files chr-by-chr: missing variants and pair LoF. R loads the results."""
import argparse
import csv
import os
import sys

def iter_chroms(plan_dir):
    plans=[f for f in os.listdir(plan_dir) if f.endswith('.mnp_gt.plan.txt')]
    chroms=sorted(set(f.replace('chr','').replace('.mnp_gt.plan.txt','') for f in plans))
    return chroms

def load_plan(path):
    rows=[]
    with open(path,newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            if len(row)<8:
                continue
            rows.append(row)
    return rows

def load_vep_ids_by_sample(path):
    """Return set of (sample_id, variant_id) present in VEP CSV."""
    seen=set()
    with open(path,newline='') as f:
        r=csv.DictReader(f)
        if 'ID' not in r.fieldnames:
            return seen
        for row in r:
            sid=(row.get('Sample.ID') or '').strip()
            vid=(row.get('ID') or '').strip()
            if sid and vid:
                seen.add((sid,vid))
    return seen

def load_vep_first_per_id(path,id_col='ID',lof_col='Variant.LoF_level'):
    """Return dict variant_id -> (lof_value, row_dict). First row per ID."""
    first={}
    with open(path,newline='') as f:
        r=csv.DictReader(f)
        for row in r:
            vid=(row.get(id_col) or '').strip()
            if not vid or vid in first:
                continue
            lof=(row.get(lof_col) or '').strip()
            first[vid]=(lof,row)
    return first

def main():
    ap=argparse.ArgumentParser(description='Build MNP report summary files for R (chr-by-chr).')
    ap.add_argument('-p','--plan-dir',default='data/mnp',help='Directory with chr*.mnp_gt.plan.txt')
    ap.add_argument('-c','--csv-dir',default='data/mnp',help='Directory with chr*.mnp_gt.vep.csv')
    ap.add_argument('-o','--output-dir',default='data/mnp',help='Where to write summary files')
    args=ap.parse_args()

    chroms=iter_chroms(args.plan_dir)
    if not chroms:
        print('No plan files found.',file=sys.stderr)
        sys.exit(1)

    missing_rows=[]
    pair_lof_rows=[]
    n_carrier_combos_total=0
    n_carrier_combos_with_missing=0

    for ch in chroms:
        plan_path=os.path.join(args.plan_dir,f'chr{ch}.mnp_gt.plan.txt')
        vep_path=os.path.join(args.csv_dir,f'chr{ch}.mnp_gt.vep.csv')
        if not os.path.isfile(plan_path):
            continue
        if not os.path.isfile(vep_path):
            print(f'Skip chr{ch}: no VEP CSV',file=sys.stderr)
            continue

        plan=load_plan(plan_path)
        vep_ids=load_vep_ids_by_sample(vep_path)
        vep_first=load_vep_first_per_id(vep_path)

        for row in plan:
            mnp_id=(row.get('mnp_id') or '').strip()
            id1=(row.get('id1') or '').strip()
            id2=(row.get('id2') or '').strip()
            carriers_s=(row.get('carriers') or '').strip()
            chrom=(row.get('chrom') or f'chr{ch}').strip()
            pos=(row.get('pos') or '').strip()
            if not mnp_id or not id1 or not id2:
                continue
            expected_ids={id1,id2,mnp_id}

            for sample in (s.strip() for s in carriers_s.split(',') if s.strip()):
                n_carrier_combos_total+=1
                observed={vid for (sid,vid) in vep_ids if sid==sample and vid in expected_ids}
                missing=expected_ids-observed
                n_obs=len(observed)
                n_miss=len(missing)
                if n_miss>0:
                    n_carrier_combos_with_missing+=1
                    missing_rows.append({
                        'mnp_id':mnp_id,'chrom':chrom,'sample':sample,
                        'n_observed':n_obs,'n_missing':n_miss,
                        'missing_ids_str':';'.join(sorted(missing))
                    })

            lof1=vep_first.get(id1,('',''))[0]
            lof2=vep_first.get(id2,('',''))[0]
            lofm=vep_first.get(mnp_id,('',''))[0]
            mnp_path_snps_not=(lofm=='1' and (lof1!='1' or lof2!='1'))
            snps_path_mnp_not=((lof1=='1' or lof2=='1') and lofm!='1')
            pair_lof_rows.append({
                'mnp_id':mnp_id,'chrom':chrom,'pos':pos,
                'id1':id1,'id2':id2,'id_mnp':mnp_id,
                'lof_SNP1':lof1,'lof_SNP2':lof2,'lof_MNP':lofm,
                'mnp_pathogenic_snps_not':str(mnp_path_snps_not),
                'snps_pathogenic_mnp_not':str(snps_path_mnp_not)
            })

    os.makedirs(args.output_dir,exist_ok=True)
    missing_path=os.path.join(args.output_dir,'mnp_report_missing_variants.txt')
    pair_path=os.path.join(args.output_dir,'mnp_report_pair_lof.txt')
    stats_path=os.path.join(args.output_dir,'mnp_report_summary_stats.txt')

    with open(missing_path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['mnp_id','chrom','sample','n_observed','n_missing','missing_ids_str'],delimiter='\t')
        w.writeheader()
        w.writerows(missing_rows)

    with open(pair_path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['mnp_id','chrom','pos','id1','id2','id_mnp','lof_SNP1','lof_SNP2','lof_MNP','mnp_pathogenic_snps_not','snps_pathogenic_mnp_not'],delimiter='\t')
        w.writeheader()
        w.writerows(pair_lof_rows)

    with open(stats_path,'w') as f:
        f.write(f'n_carrier_combos_total\t{n_carrier_combos_total}\n')
        f.write(f'n_carrier_combos_with_missing\t{n_carrier_combos_with_missing}\n')

    print(f'Wrote {missing_path} ({len(missing_rows)} rows with missing)',file=sys.stderr)
    print(f'Wrote {pair_path} ({len(pair_lof_rows)} pairs)',file=sys.stderr)
    print(f'Wrote {stats_path}',file=sys.stderr)

if __name__=='__main__':
    main()
