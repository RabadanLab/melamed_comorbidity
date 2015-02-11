f = open('census/census_germline')
cens = f.read().split("\n")
cens_germline = list(set(cens) & set(symbol_info.keys()))

miss = set(cens) - set(symbol_info.keys())

for m in miss:
    if m in alias_to_symbol:
        if len(alias_to_symbol[m]) == 1:
            cens_germline.append(alias_to_symbol[m][0])
        else:
            print m + ' has > 1'
    else:
        print m + ' is not there'

# manual review...
cens_germline.append('FLCN')
f = open('census/germline_to_remove','w')
f.write('\n'.join(cens_germline)+ '\n')
f.close()
