#!/usr/bin/env python
"""
Running `elpmpp02.py main` generates `ELPMPP02.PY.TXT`. This script transforms the output file 
to test cases.
"""
import os

PATH_SRC = os.path.realpath(os.path.dirname(__file__))

def main():
    with open(os.path.join(PATH_SRC, 'ELPMPP02.PY.TXT')) as f:
        lines = f.readlines()

    cases: dict[float, list] = {}

    for line in lines:
        if line.find('JD') > -1:
            cur_t = float(line.strip().split()[-1])
            if cur_t not in cases:
                cases[cur_t] = []
            continue
        elif len(line.strip()) > 0:
            is_fp = False
            try:
                to_remove = ['X', 'Y', 'Z', 'km/day', 'km', '=', '\'']
                for s in to_remove:
                    line = line.replace(s, '')
                line = line.strip()
                values = [float(e) for e in line.split()]
                is_fp = len(values) == 3
            except ValueError:
                values = None
            if is_fp:
                cases[cur_t].extend(line.split())

    case_code = ''

    dates = sorted(cases.keys())
    case_code +=f'    #[test]\n'
    case_code +=f'    fn check_values() {{\n'
    case_code += '        let cases = vec![\n'
    for cur_date in dates:
        if cur_date < 2400000.0 or cur_date > 2460000.0:
            continue
        case_code += f'            ({cur_date:.1f}, \n'
        case_code += '                (' + ', '.join(cases[cur_date]) + '),\n'
        case_code += '            ),\n'
    case_code += '        ];\n'
    case_code += f"""
        let mut failcnt = 0;
        for (t, expected) in cases.iter() {{
            let actual = cartesian(*t, 0.0);

            let mut failed = false;
            if !is_matched(expected, &actual) {{
                failed = true;
            }}
            if !failed {{
                println!("[t={{t}}] PASSED")
            }} else {{
                failcnt += 1; 
                println!("[t={{t}}] FAILED")
            }}
        }}

        assert_eq!(failcnt, 0);\n"""
    case_code += '    }\n'

    print(case_code)

if __name__ == '__main__':
    main()