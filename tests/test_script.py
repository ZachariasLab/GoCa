import io
import os

os.chdir('tests')

def test_1a0r_1():
    os.system('echo "2\n3\n" | ../GoCa 1a0r_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('1a0r_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('1a0r_test1.gro'))

def test_1uao_1():
    os.system('../GoCa 1uao_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('1uao_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('1uao_test1.gro'))

def test_1uao_2():
    os.system('../GoCa 1uao_test2.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('1uao_test2.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('1uao_test2.gro'))

def test_1uao_3():
    os.system('../GoCa 1uao_test3.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('1uao_test3.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('1uao_test3.gro'))

def test_2fka_1():
    os.system('../GoCa 2fka_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('2fka_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('2fka_test1.gro'))

def test_2ma2_1():
    os.system('../GoCa 2ma2_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('2ma2_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('2ma2_test1.gro'))
    assert list(io.open('tables/table_d0.xvg')) == list(io.open('2ma2_test1_table_d0.xvg'))

def test_3chy_1():
    os.system('../GoCa 3chy_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('3chy_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('3chy_test1.gro'))

def test_6ezj_1():
    os.system('../GoCa 6ezj_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('6ezj_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('6ezj_test1.gro'))

def test_6bno_1():
    os.system('echo "2 3\n" | ../GoCa 6bno_test1.ini')
    assert list(io.open('output.top'))[2:] == list(io.open('6bno_test1.top'))[2:]
    assert list(io.open('output.gro')) == list(io.open('6bno_test1.gro'))

