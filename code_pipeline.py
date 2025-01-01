
from collections import defaultdict
import csv

# Fonction pour lire un fichier FASTA et retourner une liste de séquences
def read_fasta(filename):
    sequences = []
    with open(filename, 'r') as f:
        sequence = ""
        for line in f:
            if line.startswith(">"):  # Ignore les titres
                if sequence:
                    sequences.append(sequence)
                sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

# Fonction pour extraire les k-mers d'une séquence donnée
def extract_kmers(sequence, k=21):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

# Fonction pour créer un dictionnaire de k-mers et leur fréquence
def count_kmers(sequences, k=21):
    kmer_dict = defaultdict(int)
    for seq in sequences:
        kmers = extract_kmers(seq, k)
        for kmer in kmers:
            kmer_dict[kmer] += 1
    return kmer_dict

# Fonction pour charger les informations du fichier documentaire
def load_documentary(file_documentary):
    kmers_info = {}
    with open(file_documentary, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Sauter l'entête
        for row in reader:
            kmer, kmer_type = row[0], row[1]
            kmers_info[kmer] = kmer_type
    return kmers_info

# Fonction pour comparer les k-mers entre les séquences et l'assemblage
def compare_kmers(dic_sequence, dic_assemblage, kmers_info, output_file, score_file):
    # Calcul de la matrice et du score
    total_diff = 0
    total_seq_count = 0
    total_compression = 0
    total_ajout = 0
    
    with open(output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(["21-mer", "Type", "Compte séquences", "Compte assemblage", "Différence"])
        
        for kmer, kmer_type in kmers_info.items():
            count_seq = dic_sequence.get(kmer, 0)
            count_assemblage = dic_assemblage.get(kmer, 0)
            diff = count_seq - count_assemblage
            
            # Écrire la ligne dans le fichier comparatif
            writer.writerow([kmer, kmer_type, count_seq, count_assemblage, diff])
            
            # Calcul du score
            total_diff += abs(diff)
            total_seq_count += count_seq
            if diff < 0:
                total_compression += abs(diff)
            elif diff > 0:
                total_ajout += diff
    
    # Calcul du score global
    score = total_diff / total_seq_count if total_seq_count != 0 else 0
    compression = total_compression / total_seq_count if total_seq_count != 0 else 0
    ajout = total_ajout / total_seq_count if total_seq_count != 0 else 0
    
    with open(score_file, 'w') as score_f:
        score_f.write(f"Score: {score}\n")
        score_f.write(f"Compression: {compression}\n")
        score_f.write(f"Ajout: {ajout}\n")

# Main
if __name__ == "__main__":
    # Fichiers d'entrée
    sequences_file = "sequences.fasta"
    assemblage_file = "assemblage.fasta"
    file_documentary = "documentary.csv"

    # Sorties
    output_file = "comparatif.csv"
    score_file = "score.txt"

    # Lecture des séquences et de l'assemblage
    sequences = read_fasta(sequences_file)
    assemblage = read_fasta(assemblage_file)

    # Comptage des k-mers dans les séquences et l'assemblage
    dic_sequence = count_kmers(sequences)
    dic_assemblage = count_kmers(assemblage)

    # Chargement des informations du fichier documentaire
    kmers_info = load_documentary(file_documentary)

    # Comparaison et calcul des scores
    compare_kmers(dic_sequence, dic_assemblage, kmers_info, output_file, score_file)

