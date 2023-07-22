cd bed
ls *.bed |
while read file_name;
do
cat "$file_name" >> all.bed
echo "\n" >> all.bed
done