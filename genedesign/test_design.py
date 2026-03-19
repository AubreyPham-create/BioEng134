import os
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.hairpin_checker import hairpin_checker

def verify_final_transcript(designer, transcript):
    """
    Step 2: Takes a completed Transcript object and runs it through the checkers.
    """
    # Extract data from the transcript object produced by the designer
    gene_dna = "".join(transcript.codons).upper()
    peptide_target = transcript.peptide
    
    print(f"\n--- Step 2: Auditing DNA for {len(peptide_target)} AA Protein ---")
    print("="*30)
    print("   GENE BODY AUDIT (Excluding RBS)")
    print("="*30)
    print(f"Gene Length: {len(gene_dna)} bp")
    
    # 1. Check Forbidden Sequences (Restriction Sites)
    passed_forbidden, forbidden_site= designer.forbidden_checker.run(gene_dna)
    print(f"Forbidden Sequence Free: {'✅' if passed_forbidden else '❌'}")

    # 2. Check Internal Promoters
    passed_promoter, found_promoter = designer.promoter_checker.run(gene_dna)
    print(f"Internal Promoter Free:  {'✅' if passed_promoter else '❌'}")

    # 3. Check Homopolymers
    is_homopolymer_free = designer.homopolymer_checker.run(gene_dna)
    print(f"Homopolymer Free:        {'✅' if is_homopolymer_free else '❌'}")

    # 4. Check for Hairpins (Secondary Structure)
    passed_hairpin, hairpin_info = hairpin_checker(gene_dna)
    print(f"Hairpin Free (Stable):   {'✅' if passed_hairpin else '❌'}")
    
    # 5. Translation Match (Verification)
    codon_to_aa = {}
    for aa, data in designer.codon_data.items():
        for codon, weight in data:
            codon_to_aa[codon] = aa
            
    translated = "".join([codon_to_aa.get(c, "?") for c in transcript.codons])
    
    # Use rstrip('*') to ensure comparison works even if stop codons are added
    if translated.rstrip('*') == peptide_target.rstrip('*'):
        print(f"Translation Match:       ✅")
    else:
        print(f"Translation Match:       ❌ (Mismatch detected!)")

def test_designer():
    designer = TranscriptDesigner()
    designer.initiate()
    
    test_peptide = (
        "MKIVTYSKNVFIPLTRNCRNRCDYCGFRSEEIGVMSAEEVRKILSAARGAKEALFTFGERPDEVYPEFKAMLLE"
        "MGYKSMVDYILEMNKIAVEMGFLPHTNAGILSKEEMRKLKPYNASMGLMLEQAVELECHANSPGKKPEVRIKM"
        "IRDAGKLQIPFTTGILVGIGEGREDRLYSLEIIAEIQDNYGHIQEVIVQNFKPKKGTPMENTPPPTLEEMLEAV"
        "KAAREILPQEVAIQIPPNLVDDIYPFLKAGANDVGGISNVTHDYINPESPWPEVERLERALRGEFILKERLPIY"
        "PRFVKLKWYGEALEPLIEKYSDADGYARP"
        "MRVEVFPVEGLPLIKEGDDLAELISSRVRFEDGDVLVVCSTVISKAEGRIRRLEEFNPSE"
        "RAKEIAARIGKPAEFVQAVLEESEEVLLDFPFLLVKAKFGNVCVNAGIDASNVEEGSLLL"
        "PPLDPDGSAEKLRRRILELTGKRVGVIITDTNGRCFRRGVVGFAIGISGVKAMKDWIGRK"
        "DLYGRELEVTVECVADEIAAFANLLMGEGGDGIPAVVVRGLNVAGEGSMEEIYRSEEEDV"
        "IRRCLKRCL"
    )
    
    print(f"--- Starting Design for {len(test_peptide)} AA Protein ---")
    transcript = designer.run(test_peptide, ignores=set())
    verify_final_transcript(designer, transcript)

def test_new_peptide():
    designer = TranscriptDesigner()
    designer.initiate()
    
    new_peptide = (
        "MEVAFDIRNFSVYYGNKVGIRNVNLEIYRNKVTAIIGPSGCGKSTLLRSLNRLIELVDGV"
        "RIEGEVIFDGKNIYDDDVDAVELRRRIGMVFQHPNPFPKSIFDNVAYGPRVHGIKDKERL"
        "KEIVEESLKKAALWDEVKDRLNDSALGLSGGQQQRLCIARAIATNPEVILFDEPTSALDP"
        "IAAAKIEELMVELKKNYTVVVVTHNIQQAARISDYVAFFWMGELVEYGKTAKVFEKPEKE"
        "LTEKYLTGRVG"
    )
    
    print(f"--- Running Design for New Peptide ({len(new_peptide)} AA) ---")
    transcript = designer.run(new_peptide, ignores=set())
    verify_final_transcript(designer, transcript)

def test_membrane_protein():
    designer = TranscriptDesigner()
    designer.initiate()
    
    membrane_peptide = (
        "MLGSLVMIGVGLLLVYLGIVKKMEPLLLVPIGIGAILVNIPGGGLAEEGSIFDLFLKYLI"
        "HTEIVPLLIFLGLGALTDFSPLLANPKTFLLGAAAQIGIFAALIAALFLGFTPQEAASIG"
        "IIGGADGPTTIYTTTILAPHLLAATAVAAYSYMSLVPIIQPPIIKALTSSRERKIKMRQL"
        "RIVSKKEKILFPIATIIISGFLAPKALPLVGMLMTGNLFRESGVTDRLAKGASEELMNIM"
        "TIILGLSVGSTMRAESFLTQKTLLVLALGVVAFAAATAGGVLLAKVMNLFLKEKINPMIG"
        "AAGVSAVPMSARVVQRLAIEEDPHNHILMHAMGPNVAGVIGSAVAAGVLIQILG"
    )
    
    print(f"--- Designing 900bp Membrane Protein Sequence ---")
    transcript = designer.run(membrane_peptide, ignores=set())
    verify_final_transcript(designer, transcript)

def test_large_enzyme():
    # 1. Initialize
    designer = TranscriptDesigner()
    designer.initiate()
    
    # 2. 660 AA Peptide Sequence
    large_peptide = (
        "MYSVKWDNKTGGIKLLPSGDGDVASMVRPVFYEELDILGFGKYWDYPKSEEPLLWSVGRS"
        "YFYNGEKVAQASGGGFYEYPKLHIFKEGLSIEPVDVDKMVKENIEIIESVAHTTLDFISE"
        "MYVKYKNKIDKTVVAFSGGKDSTVLLDLVQRVLSPDEYIVVFNDTTMELSPTYEFIEKIK"
        "KKYSNLMFLVSRYEKPAIEMWREVGVPSKIHRWCCTIYKTVPTIKLIQELIKKETPKILL"
        "YDGIRGDESLTRAKLSHISRGKHLQQINIHPIIDWNSAMVYLYIFMRDLPLNKLYRYGIV"
        "RVGCAVCPFESYWKEAILWLKFRDEIKPYLEVIEKYARDRRCKGSELKEFVKEGVWKKRV"
        "GNKGLNLADKVNITEEKDKATVISINTPNENWFEWIKTLGKVVVGNPSGYIEIKGIQYPF"
        "KWQKTGNNIKIKFTSKIKNDDLFYIKNVSYKSAYCIRCGTCEIECPSSAIKIERDKIRID"
        "EKKCKNCYRCLTFIDRGCIVADSYRESKVVKSMNNIARYRTFGMRKEWLEEFFSDTDKWW"
        "IENKLGPIQKISMRLWLTDSEILEVKGNKIILSKTGELLREIGASDLFTWAIIWTNLVRN"
        "SPLIEWYVQRIPWGRAYNREELIQTLGEEYRVRTRKNAILELVQLFDYTPLGRELKIGVP"
        "VKRGKTTTAIEKIGLTNESDLKNLHQLAVLYSLYRYSEKTSRYHFTVSELYEASAIEGPY"
        "KLFGIPKTLLEKILIGLEELYGKNWINVELVADLDNINLNPTKNSMDVLKLYAGVEYD"
    )
    
    print(f"--- Starting Design for Large Enzyme ({len(large_peptide)} AA) ---")
    
    # 3. Run Design
    transcript = designer.run(large_peptide, ignores=set())
    
    # 4. Run Audit
    verify_final_transcript(designer, transcript)

if __name__ == "__main__":
    # Executing the full test suite
    test_designer()
    test_new_peptide()
    test_membrane_protein()
    test_large_enzyme()