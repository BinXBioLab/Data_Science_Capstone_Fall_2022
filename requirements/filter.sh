touch requirements-colab.txt
cp requirements.txt requirements-colab.txt

for i in $(cat colab-exclude.txt)
do
	sed -i "" "/^$i/d" requirements-colab.txt
done
