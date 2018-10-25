# coding=UTF-8
#定义一个数组，用于实现哈希表
class Array(object): 
    def __init__(self,size=32,init=None):
        self._size = size
        self._items = [init] * size
    def __getitem__(self, index):
        return self._items[index]
    def __setitem__(self, index, value):
        self._items[index] = value
    def __len__(self):
        return self._size
    def clear(self,value=None):
        for i in range(len(self._items)):
            self._items[i] = value
    def __iter__(self):
        for item in self._items:
            yield item

'''
定义一个hash表数组的槽
注意，一个槽有三种状态。
1.从未使用HashTable.UNUSED。此槽没有被使用和冲突过，查找时只要找到UNUSED就不用再继续探查了
2.使用过但是remove了，此时是HashTable.EMPTY，该探查点后边的元素扔可能是有key
3.槽正在使用Slot节点
'''
class Slot(object):
    def __init__(self,key,value):
        self.key = key
        self.value = value
#定义一个哈希表类
class HashTable(object):
    UNUSED = None #槽没被使用过(槽指的是数组中的一个位置)
    EMPTY = Slot(None,None)#使用过却被删除了

    def __init__(self):
        self._table = Array(8,init=HashTable.UNUSED)#初始化槽没被使用过
        self.length = 0
    @property
    def _load_factor(self):#装载因子(load factor)，当装载因子超过某一阈值时，需要进行重哈希操作
        return self.length/float(len(self._table))
    def __len__(self):#哈希表的长度
        return self.length
    def _hash(self,key):#定义hash函数，返回哈希表数组的索引
        return abs(hash(key))%len(self._table)
    def _find_key(self,key):#搜寻key键，如果存在，返回key键对应的索引，否则返回None
        index = self._hash(key)
        _len = len(self._table)
        while self._table[index] is not HashTable.UNUSED:
            if self._table[index] is HashTable.EMPTY:
                index = (index*5 +1)%_len#解决哈希冲突的一种方式
                continue
            elif self._table[index].key == key:
                return index
            else:
                index = (index*5 + 1)%_len
        return None
    def _find_slot_for_insert(self,key):#寻找哈希表中可以插入的槽，
        index = self._hash(key)
        _len = len(self._table)
        while not self._slot_can_insert(index):
            index = (index * 5 + 1)%_len#解决哈希冲突的一种方式
        return  index
    def _slot_can_insert(self,index):#判断发现的槽能否插入，空槽或从未被使用即返回True
        return (self._table[index] is HashTable.EMPTY or self._table[index] is HashTable.UNUSED)
    def __contains__(self, key):#in 操作符，判断key是否在hash表里，如果存在返回其索引，否则返回None
        index = self._find_key(key)
        return index is not None
    def add(self,key,value):#添加元素操作，有 key 则更新，否则插入
        if key in self:#调用__contains__判断key是否在哈希表中
            index = self._find_key(key)#搜寻key的索引
            self._table[index].value = value#更新哈希表中key对应的value
            return False
        else:#如果key不再hash表中，则添加元素Slot(key,value),并更新长度
            index = self._find_slot_for_insert(key)
            self._table[index] = Slot(key,value)
            self.length += 1
            if self._load_factor >= 0.8:#如果装载因子大于0.8，则进行重哈希操作
                self._rehash()
            return True

# '''
# 重哈希(Rehashing)，
# 步骤就是重新开辟一块新的空间。本函数里的策略是扩大为已经使用的槽数目的两倍。
# 开辟了新空间以后，会把原来哈希表里不为空槽的数据重新插入到新的哈希表里，
# 插入方式和之前一样。这就是 rehashing 操作。
# '''
    def _rehash(self):
        old_table = self._table
        newsize = len(self._table) * 2
        self._table = Array(newsize,HashTable.UNUSED)#定义新的哈希表
        self.length = 0#初始长度为0
        for slot in old_table:#把原来哈希表里不为空槽的数据重新插入到新的哈希表里
            if slot is not HashTable.UNUSED and slot is not HashTable.EMPTY:
                index = self._find_slot_for_insert(slot.key)
                self._table[index]=slot
                self.length += 1
    def get(self,key,default = None):#获取 key 的值，不存在返回默认值 None
        index = self._find_key(key)
        if index is None:
            return default
        else:
            return self._table[index].value
    def remove(self,key):#删除一个 key，这里其实不是真删除，而是标记为 Empty
        index = self._find_key(key)
        if index is None:
            raise KeyError()
        value = self._table[index].value
        self.length -= 1
        self._table[index] = HashTable.EMPTY
        return value
    def __iter__(self):
        for slot in self._table:
            if slot not in (HashTable.EMPTY,HashTable.UNUSED):
                yield slot.key