# 绘制每种可能交集的情况，左图为类别元素个数，下图为所有类别交集种类，上图为对于种类元素个数.Cn1+Cn2+Cn3……Cnn
mammals = ["Cat", "Dog", "Horse", "Sheep", "Pig", "Cattle", "Rhinoceros", "Moose"]
herbivores = ["Horse", "Sheep", "Cattle", "Moose", "Rhinoceros"]
domesticated = ["Dog", "Chicken", "Horse", "Sheep", "Pig", "Cattle", "Duck"]
(mammals, herbivores, domesticated)
from Python.basic.upsetplot import from_contents,plot,UpSet

animals = from_contents(
    {"mammal": mammals, "herbivore": herbivores, "domesticated": domesticated}
)
animals
plot(animals)
# ax_dict = UpSet(animals, subset_size="count").plot()