for i in "2.5,2" "2.5,2.5"
do
echo $i
python build_network.py $i
python run_save_network.py $i
done
