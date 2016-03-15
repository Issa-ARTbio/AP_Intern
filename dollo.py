def dollo(node):
    """Infer the state of domains for all ancestral nodes in a given tree,
    under the assumption that a domain gain is so expensive that it can only be
    gained once.  Also calculate the gain and loss events for all nodes on the
    tree.
    This function expects a bifurcating tree as an ete2-object, in which the
    leafs contain the feature 'domains' (a set of all domains present in the
    leaf).
    """
    leaves_to_root(node)
    root_to_leaves(node)
    events(node)


def leaves_to_root(node):
    """For one node, collect all the domains that are present in the children.
    If the domains are present in all children, set the state for this domain
    to 1, if a domain is only present in some of the children, set the state to
    0.  After running this, the root node will have a status for all domains in
    the tree.
    """
    if node.is_leaf():
        return dict((domain, 1) for domain in node.domains)

    state = dict()
    call = node.get_children()
    lstate = [ ]
    for child in call :
        lstate.append( leaves_to_root( child ) )
    ldomains = [ set( s.keys() ) for s in lstate ]
    dom_inter, dom_union = set(), set()
    if len(ldomains) > 0:
        dom_inter = ldomains[0]
        dom_union = ldomains[0]
    for sdomain in ldomains[1:]:
        dom_inter = dom_inter.intersection(sdomain)
        dom_union = dom_union.union(sdomain)
    for domain in dom_union :
        if domain in dom_inter:
            state[domain] = 1  # this domain is definitely there
        else:
            state[domain] = 0  # domain status is unknown at this node

    node.add_features(domains=state)
    return node.domains


def root_to_leaves(node):
    """Determine the state for uncertain domains in a node according to the
    state it has in the parents.  If a domain status was 0 in a node, but is
    available in the parent, it must have been present in the parent (a domain
    can only be gained once).
    """
    if node.is_leaf():
        return
    elif node.is_root():
        # delete all domains with count 0 at root and turn the dict into a set
        # == keep all domains with count 1
        node.domains = set(d for (d, c) in node.domains.items() if c == 1)
    else:
        # delete all domains with count 0 if the domain does not exist in the
        # parent.
        item_node = list(node.domains.items())
        for d, c in item_node:  # work on a copy
            if c == 0 and d not in node.up.domains:
                del node.domains[d]
        # count numbers are now irrelevant
        node.domains = set(node.domains)
    for child in node.get_children( )  :
        root_to_leaves( child )

def events(node):
    """Calculate which domains have been gained or lost on all nodes in the
    tree.
    """
    d = node.domains
    if node.is_root():
        p = set()
    else:
        p = node.up.domains
    gained = d.difference(p)
    lost = p.difference(d)
    node.add_features(gained_domains=gained, lost_domains=lost)

    if node.is_leaf():
        return

    for child in node.get_children( ) :
        events( child )

if __name__ == '__main__':
    import ete3
    # import dollo
    import sys
    # ete3 can also directly read a file
    tree = ete3.Tree(sys.argv[1])
#     tree = ete3.Tree('(((Synechococcus_sp_PCC_7336:0.31351706318466310286,(Cyanobacterium_YellowstoneA:0.05103864297600101824,Cyanobacterium_YellowstoneB:0.04567157612933206434):0.26732139512279695648):0.11651441969273111654,((((Synechococcus_sp_PCC_7502:0.22557752665191099783,Pseudanabaena_sp_PCC_7429:0.24288852599962218459):0.03473507036651951596,Pseudanabaena_sp_PCC_6802:0.17858438820673211422):0.04214269083546344496,Pseudanabaena_sp_PCC_7367:0.24435547467613397132):0.16511126852399632403,(Gloeomargarita_lithophora:0.52273917959513516163,((((Synechococcus_sp_PCC_6312:0.21236223177778931759,(Thermosynechococcus_elongatus_BP1:0.21423572777315644244,Synechococcus_calcipolaris:0.16965407183620034859):0.06425385377127000586):0.10094610025167553846,Cyanothece_sp_PCC_7425:0.18502732626359114088):0.04490744762183230404,(Acaryochloris_marina_MBIC11017:0.01217969585970165095,Acaryochloris_sp_CCMEE_5410:0.01822215686830522957):0.27426610433018555613):0.04545653694923947052,((((Synechococcus_elongatus_PCC7942:0.00005313973350376229,Synechococcus_elongatus_PCC6301:0.00230870286914624861):0.23097894378750524758,(Synechococcus_sp_RCC307:0.20976955763602134208,(Synechococcus_sp_WH5701:0.20576440303378293328,(((Synechococcus_sp_CC9311:0.14029393788446364866,Synechococcus_sp_WH7803:0.11800010033254393349):0.05587345873787614992,(Synechococcus_sp_WH8102:0.09632283481312754747,(Synechococcus_sp_CC9605:0.09523614251108808437,Synechococcus_sp_CC9902:0.12150486637192933759):0.02515429990853526934):0.07711151291003463804):0.03048673168143171616,((Prochlorococcus_marinus_CCMP1375:0.24669308921211188790,((Prochlorococcus_marinus_NATL1A:0.01034134087011697248,Prochlorococcus_marinus_NATL2A:0.01039148682608429138):0.23535849937188343950,((Prochlorococcus_marinus_MIT9515:0.06023895379249264576,Prochlorococcus_marinus_CCMP1986:0.04854230525455979078):0.06372866097907793625,(Prochlorococcus_marinus_MIT9312:0.03360214040564027393,((Prochlorococcus_marinus_MIT9301:0.02671299519911931805,Prochlorococcus_marinus_AS9601:0.02289906692463827831):0.00770378656974180995,Prochlorococcus_marinus_MIT9215:0.03928354067506659858):0.01622692955547334198):0.07209062957838605068):0.42066605440976628794):0.07903894487403052838):0.17401317914634775730,(Prochlorococcus_marinus_MIT9303:0.01067854629499025154,Prochlorococcus_marinus_MIT9313:0.01269588385027749060):0.13016618253388351212):0.06961724052428064358):0.09677713345840242842):0.05026578438065273935):0.47715510146078882192):0.12923105205157645048,((Synechococcus_sp_PCC_7335:0.24712603382396730600,Leptolyngbya_sp_PCC_7375:0.20511757464869290191):0.08144255874183710386,(Leptolyngbya_sp_PCC_6406:0.19346310961520099547,Nodosilinea_nodulosa_PCC_7104:0.20958223568088149569):0.04844921642989390848):0.06935263480185017981):0.03234659282676372732,(((((Rubidibacter_lacunae_KORDI_512:0.28698483477864222824,(Halothece_sp_PCC_7418:0.07272029524926727773,Dactylococcopsis_salina_PCC_8305:0.08824397904242899104):0.17061147510513621772):0.07834384490573383097,((Spirulina_subsalsa_PCC_9445:0.16537184694608739188,Spirulina_major_PCC_6313:0.22997714169734409517):0.08783168623251645657,((Gloeocapsa_sp_PCC_73106:0.27442727072068212602,((Leptolyngbya_sp_PCC_7376:0.12268257252114714295,Synechococcus_sp_PCC_7002:0.09612027801542420702):0.18842895248434690658,(Geminocystis_herdmanii_PCC_6308:0.14432757194691950287,Cyanobacterium_stanieri_PCC_7202:0.14362431064063663211):0.15492081715972808031):0.03747079698934641101):0.02075837059739236642,((Stanieria_sp_PCC_7437:0.12318946826573846931,(Pleurocapsa_sp_PCC_7319:0.16566350568379092922,Chroococcidiopsis_sp_PCC_6712:0.16310117114117184123):0.03855557815779731695):0.07703685776432613042,(Pleurocapsa_sp_PCC_7327:0.15119197086110203188,((Cyanothece_sp_PCC_8801:0.11528200432560758992,(((Cyanothece_sp_ATCC51142:0.00013213943730671607,Cyanothece_sp_ATCC_51472:0.00010062599064728421):0.02857537512636421431,Cyanothece_sp_CCY0110:0.03420861716435349281):0.02539245395556166338,Crocosphaera_watsonii_WH8501:0.07773142363495394447):0.11097394450644292030):0.06306702913632684926,(Synechocystis_sp_PCC6803:0.
# 28377217331711956927,(Microcystis_aeruginosa_PCC_7806:0.01466399804614199401,Microcystis_aeruginosa_NIES843:0.00913370343134168802):0.19937501093122780849):0.03476403102142960194):0.02105882814723538179):0.02986432097994805301):0.01721903061672224941):0.02985674609995001005):0.02106499837831399269):0.04619965014263693881,(Moorea_producens_3L:0.17196742120983363189,(Coleofasciculus_chthonoplastes_PCC_7420:0.15174853361990467415,Microcoleus_sp_PCC_7113:0.12185129820001944223):0.02324787073661519279):0.04561585872812332770):0.03075052831930170927,(((Chroococcidiopsis_thermalis_PCC_7203:0.15626904723940410191,(Synechocystis_sp_PCC_7509:0.18269106917880961749,Gloeocapsa_sp_PCC_7428:0.12029691284599618173):0.02075077178344942111):0.02938118188613386347,((Tolypothrix_sp_PCC_9009:0.10617039295284882994,((Nodularia_spumigena_CCY9414:0.09962480895188752239,(Nostoc_sp_PCC_7107:0.08112012777183683077,(Nostoc_sp_PCC_7524:0.05398244441490487877,(Nostoc_sp_PCC7120:0.01206848985381956188,Anabaena_variabilis:0.01058748144122514945):0.05370655720022284724):0.02414139663273106098):0.01566242002257948671):0.01287964152782410977,((Calothrix_sp_PCC_7507:0.04412796476689943154,Microchaete_sp_PCC_7126:0.05972288102666735432):0.04165324958875344513,(Nostoc_punctiforme:0.07744440457305856729,(Cylindrospermum_stagnale_PCC_7417:0.06268127120703663457,(((Raphidiopsis_brookii_D9:0.02022900316270691973,Cylindrospermopsis_raciborskii_CS505:0.02274833233737944341):0.15275742395010846741,Trichormus_azollae_0708:0.05912398493651230269):0.01750481648985037844,(Anabaena_cylindrica_PCC_7122:0.04364492569714339965,Anabaena_sp_PCC_7108:0.04938240809882146060):0.01375885607114792596):0.04379659835895828579):0.02169818388444162693):0.00928289070973055061):0.00893655556779493140):0.02646300335615132601):0.02053568435636278153,((Richelia_intracellularis_HH01:0.27098549262331739218,(Rivularia_sp_PCC_7116:0.16522285503520173222,Calothrix_sp_PCC_6303:0.16120332687796173898):0.02022431975489523068):0.01934526394497537022,((Mastigocladopsis_repens_PCC_10914:0.06998555657456660695,Scytonema_hofmanni_PCC_7110:0.09863929578204404247):0.03161973675362823172,((cyanobacterium_PCC_7702:0.09246299678908075081,(Chlorogloeopsis_fritschii_PCC_6912:0.00106240007642685553,Chlorogloeopsis_fritschii_PCC_9212:0.00080825301927448130):0.03972852852036402177):0.02728434805394009191,(Fischerella_sp_PCC_9605:0.04478241760007761058,((Fischerella_sp_PCC_9431:0.01766349602578041675,(Fischerella_muscicola_PCC_73103:0.01382349640307813966,Fischerella_sp_PCC_9339:0.02183227387052787474):0.00443328940458307102):0.02555025221772255078,(Fischerella_muscicola_PCC_7414:0.01311248567771962012,(Fischerella_sp_JSC11:0.00029208693032646706,Fischerella_thermalis_PCC_7521:0.00037630815328235025):0.01596825561328989521):0.01562808578902291287):0.03455921126909176216):0.01998226544290006984):0.03467573777948307345):0.01268110465243722611):0.00963846778776244681):0.07505157192120910647):0.05423000109106278915,((Oscillatoria_sp_PCC_10802:0.16032491376933960536,Oscillatoria_acuminata_PCC_6304:0.20846838070526313436):0.03261980872422277911,(((Microcoleus_vaginatus_FGP2:0.01896640582407359971,Oscillatoria_sp_PCC_7112:0.01718413337967054574):0.11127214066791790703,(Oscillatoria_sp_PCC_6506:0.00094403297350733315,Oscillatoria_sp_PCC_6407:0.00016081697036633870):0.10817577097434806499):0.06975736421329398196,((Lyngbya_sp_PCC8106:0.14551556779187130308,((Arthrospira_platensis_NIES39:0.00228300544197777323,Arthrospira_platensis_Paraca:0.00128880185796440287):0.01312680210602682265,Arthrospira_maxima_CS328:0.01870198978334828427):0.16222044582450734218):0.06042390387917623285,Trichodesmium_erythraeum:0.23233572546585579133):0.03891122916771631857):0.03215566487587199346):0.02522911994984410738):0.01895056062215981843):0.03380129887396540744,(Geitlerinema_sp_PCC_7407:0.19180575818364700713,Leptolyngbya_boryana_PCC_6306:0.23233607324975219521):0.03389427103008481323):0.02618121970183569583):0.02720711992588421835):0.03882643881640874400):0.03863333779230554432):0.09280803193775634397):0.21763435107663511858,(Gloeobacter_
# violaceus_PCC7421:0.14659015636126423998,Gloeobacter_kilaueensis_JS1:0.14280631503665322524):0.21763435107663511858);')
    # tree = ete3.Tree('(((sp1,sp2),sp3),sp4);')

    # set a new attribute to each leaf of the tree
    # containing the set of domain in the corresponding species
    # (tree&'Synechococcus_sp_PCC_6312').add_features(domains=set(['A', 'B']))
    # (tree&'Synechococcus_calcipolaris').add_features(domains=set(['A', 'B']))
    # (tree&'Thermosynechococcus_elongatus_BP1').add_features(domains=set(['biom_pole', 'ATPase']))
    #
    # (tree&'Gloeomargarita_lithophora').add_features(domains=set(['biom_cyto']))
    # (tree&'Cyanothece_sp_PCC_7425').add_features(domains=set(['biom_cyto']))
    # (tree&'Chroococcidiopsis_thermalis_PCC_7203').add_features(domains=set(['biom_cyto']))

    (tree&'Synechococcus_sp_PCC_7336').add_features(domains=set(['A','B']))
    # (tree&'sp1').add_features(domains=set(['A','B']))
    # (tree&'sp2').add_features(domains=set(['B','C']))
    # (tree&'sp3').add_features(domains=set(['D','E']))
    # (tree&'sp4').add_features(domains=set(['A','E']))


    # compute dollo parsimony
    dollo(tree)

    print(tree)
    print("\ntraverse root to leaves and right to left\n")
    print("name\tdomains\tgained\tlost")
    for node in tree.traverse( ) :
        name = node.name
        if node.is_root() :
            name = "root"
        print("{}\t{}\t{}\t{}".format(name, ", ".join(node.domains), ", ".join(node.gained_domains), ", ".join(node.lost_domains)))
