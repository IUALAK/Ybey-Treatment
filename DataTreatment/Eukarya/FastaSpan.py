import os

final_text = ''
for i in os.listdir('Eukarya_samples/Results'):
    print(i)
    with open(f'Eukarya_samples/Results/{i}', 'r') as text:
        final_text = final_text + text.read()

with open('final.fasta', 'w') as f:
    f.write(final_text)
    f.close()


