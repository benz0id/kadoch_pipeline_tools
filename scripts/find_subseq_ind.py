









seq = "MAVRKKDGGPNVKYYEAADTVTQFDNVRLWLGKNYKKYIQAEPPTNKSLSSLVVQLLQFQEEVFGKHVSNAPLTKLPIKCFLDFKAGGSLCHILAAAYKFKSDQGWRRYDFQNPSRMDRNVEMFMTIEKSLVQNNCLSRPNIFLCPEIEPKLLGKLKDIIKRHQGTVTEDKNNASHVVYPVPGNLEEEEWVRPVMKRDKQVLLHWGYYPDSYDTWIPASEIEASVEDAPTPEKPRKVHAKWILDTDTFNEWMNEEDYEVNDDKNPVSRRKKISAKTLTDEVNSPDSDRRDKKGGNYKKRKRSPSPSPTPEAKKKNAKKGPSTPYTKSKRGHREEEQEDLTKDMDEPSPVPNVEEVTLPKTVNTKKDSESAPVKGGTMTDLDEQEDESMETTGKDEDENSTGNKGEQTKNPDLHEDNVTEQTHHIIIPSYAAWFDYNSVHAIERRALPEFFNGKNKSKTPEIYLAYRNFMIDTYRLNPQEYLTSTACRRNLAGDVCAIMRVHAFLEQWGLINYQVDAESRPTPMGPPPTSHFHVLADTPSGLVPLQPKTPQQTSASQQMLNFPDKGKEKPTDMQNFGLRTDMYTKKNVPSKSKAAASATREWTEQETLLLLEALEMYKDDWNKVSEHVGSRTQDECILHFLRLPIEDPYLEDSEASLGPLAYQPIPFSQSGNPVMSTVAFLASVVDPRVASAAAKSALEEFSKMKEEVPTALVEAHVRKVEEAAKVTGKADPAFGLESSGIAGTTSDEPERIEESGNDEARVEGQATDEKKEPKEPREGGGAIEEEAKEKTSEAPKKDEEKGKEGDSEKESEKSDGDPIVDPEKEKEPKEGQEEVLKEVVESEGERKTKVERDIGEGNLSTAAAAALAAAAVKAKHLAAVEERKIKSLVALLVETQMKKLEIKLRHFEELETIMDREREALEYQRQQLLADRQAFHMEQLKYAEMRARQQHFQQMHQQQQQPPPALPPGSQPIPPTGAAGPPAVHGLAVAPASVVPAPAGSGAPPGSLGPSEQIGQAGSTAGPQQQQPAGAPQPGAVPPGVPPPGPHGPSPFPNQQTPPSMMPGAVPGSGHPGVAGNAPLGLPFGMPPPPPPPAPSIIPFGSLADSISINLPAPPNLHGHHHHLPFAPGTLPPPNLPVSMANPLHPNLPATTTMPSSLPLGPGLGSAAAQSPAIVAAVQGNLLPSASPLPDPGTPLPPDPTAPSPGTVTPVPPPQ"

sub_seqs = {
    "R499": "LAGDVCAIMR",
}

poss = []
for org_pos in sub_seqs:
    print(org_pos, ':', sub_seqs[org_pos])
    try:
        ind = seq.index(sub_seqs[org_pos])
        ind += len(sub_seqs[org_pos]) - 1
    except ValueError:
        print(org_pos + ' could not be found.')
    else:
        print(org_pos + ' -> ' + seq[ind] + str(ind + 1))
        poss.append(ind)
print(','.join([str(p + 1) for p in poss]))






